#!/usr/bin/env python3
# coding:utf-8

"""
PANORAMA clustering module for pangenome gene family clustering.

This module provides functionality to cluster gene families across multiple pangenomes
using MMseqs2 with support for both linclust (fast) and cluster (sensitive) methods.
It handles sequence processing, clustering execution, and result formatting with
comprehensive error handling and progress tracking.
"""

# default libraries
from __future__ import annotations
import argparse
import logging
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, Union, Optional, List
from multiprocessing import Manager, Lock
from time import time
from shutil import rmtree
from dataclasses import dataclass, field

# installed libraries
import pandas as pd

# local libraries
from panorama.utils import mkdir
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.format.read_binaries import load_pangenomes
from panorama.alignment.utils import (
    write_pangenomes_families_sequences,
    createdb,
    PangenomeProcessingError,
    MMSeqsError,
)


logger = logging.getLogger("PANORAMA")


class ClusteringMethod:
    """
    Defines a class representing clustering methods and their choices.
    """

    CHOICES = ["linclust", "cluster"]
    LINCLUST = "linclust"
    CLUSTER = "cluster"


# Configuration constants
@dataclass(frozen=True)
class ClusteringConfig:
    """Configuration constants for clustering operations."""

    # Column names for clustering results
    CLUSTER_COLUMN_NAMES: List[str] = field(
        default_factory=lambda: ["cluster_id", "referent", "in_clust"]
    )

    # Default clustering parameters
    DEFAULT_THREADS: int = 1
    DEFAULT_IDENTITY: float = 0.5
    DEFAULT_COVERAGE: float = 0.8
    DEFAULT_COV_MODE: int = 0
    DEFAULT_EVAL: float = 1e-3

    # MMseqs2 default parameters
    DEFAULT_MMSEQS2_OPTIONS: Dict[str, Union[int, float, str]] = field(
        default_factory=lambda: {
            "comp_bias_corr": 1,
            "kmer_per_seq": 20,
            "identity": 0.5,
            "coverage": 0.8,
            "cov_mode": 0,
            "eval": 1e-3,
            "align_mode": 3,
            "max_seq_len": 65535,
            "max_reject": 2147483647,
            "clust_mode": 0,
            # Cluster-specific options
            "max_seqs": 300,
            "min_ungapped": 30,
            "sensitivity": 5.7,
        }
    )

    # Output filename
    CLUSTERING_OUTPUT: str = "clustering_results.tsv"


# Global configuration instance
CONFIG = ClusteringConfig()


class ClusteringError(Exception):
    """Custom exception for clustering-related errors."""

    pass


class ClusteringValidationError(ClusteringError):
    """Custom exception for clustering parameter validation errors."""

    pass


def _validate_clustering_parameters(
    threads: int, mmseqs2_options: Dict[str, Union[int, float, str]]
) -> None:
    """
    Validate clustering parameters.

    Args:
        threads: Number of threads (positive integer).
        mmseqs2_options: Dictionary containing MMseqs2 parameters.

    Raises:
        ClusteringValidationError: If any parameter is invalid.
    """
    if not isinstance(threads, int) or threads <= 0:
        raise ClusteringValidationError(
            f"Threads must be positive integer, got: {threads}"
        )

    if not isinstance(mmseqs2_options, dict):
        raise ClusteringValidationError("MMseqs2 options must be a dictionary")

    # Validate required MMseqs2 parameters
    required_params = ["identity", "coverage", "cov_mode"]
    for param in required_params:
        if param not in mmseqs2_options:
            raise ClusteringValidationError(
                f"Missing required MMseqs2 parameter: {param}"
            )

    # Validate parameter ranges
    identity = mmseqs2_options.get("identity", 0)
    coverage = mmseqs2_options.get("coverage", 0)
    cov_mode = mmseqs2_options.get("cov_mode", 0)

    if not isinstance(identity, (int, float)) or not (0.0 <= identity <= 1.0):
        raise ClusteringValidationError(
            f"Identity must be between 0.0 and 1.0, got: {identity}"
        )

    if not isinstance(coverage, (int, float)) or not (0.0 <= coverage <= 1.0):
        raise ClusteringValidationError(
            f"Coverage must be between 0.0 and 1.0, got: {coverage}"
        )

    if not isinstance(cov_mode, int) or not (0 <= cov_mode <= 5):
        raise ClusteringValidationError(
            f"Coverage mode must be between 0 and 5, got: {cov_mode}"
        )


def _validate_directory_access(
    directory: Path, create_if_missing: bool = False
) -> None:
    """
    Validate directory exists and is accessible.

    Args:
        directory: Directory path to validate.
        create_if_missing: Whether to create a directory if it doesn't exist.

    Raises:
        ClusteringValidationError: If directory validation fails.
    """
    if not isinstance(directory, Path):
        raise ClusteringValidationError(f"Expected Path object, got: {type(directory)}")

    if not directory.exists():
        if create_if_missing:
            try:
                directory.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                raise ClusteringValidationError(
                    f"Cannot create directory {directory}: {e}"
                )
        else:
            raise ClusteringValidationError(f"Directory does not exist: {directory}")

    if not directory.is_dir():
        raise ClusteringValidationError(f"Path is not a directory: {directory}")


def check_cluster_parameters(args: argparse.Namespace) -> None:
    """
    Validate command line arguments for clustering operations.

    This function performs comprehensive validation of all clustering parameters
    provided via command line arguments, ensuring they meet the required
    constraints before proceeding with clustering operations.

    Args:
        args: Parsed command line arguments containing clustering parameters.
            Expected attributes: tmpdir and clustering-specific parameters.

    Raises:
        ClusteringValidationError: If any parameter validation fails.
        NotADirectoryError: If tmpdir is not a valid directory.
    """
    # Validate temporary directory
    if not hasattr(args, "tmpdir") or not isinstance(args.tmpdir, Path):
        raise ClusteringValidationError("tmpdir must be a valid Path object")

    try:
        _validate_directory_access(args.tmpdir, create_if_missing=True)
    except ClusteringValidationError as e:
        raise NotADirectoryError(
            f"The given path for temporary directory is not valid: {args.tmpdir}. "
            f"Error: {e}. Please check your path and try again."
        )

    logger.debug("Clustering parameters validated successfully")


def check_pangenome_cluster(pangenome: Pangenome) -> None:
    """
    Validate that a pangenome is ready for clustering operations.

    This function checks that the pangenome has been properly processed
    and contains the necessary data for gene family clustering operations.

    Args:
        pangenome: Pangenome object to validate.
            Must have clustered genes and associated gene family sequences.

    Raises:
        AttributeError: If pangenome is missing required data or processing steps.
    """
    if not hasattr(pangenome, "status"):
        raise AttributeError("Pangenome object missing status information")

    # Check if genes have been clustered
    genes_clustered = pangenome.status.get("genesClustered", "No")
    if genes_clustered == "No":
        raise AttributeError(
            f"Genes not clustered for pangenome {getattr(pangenome, 'name', 'unknown')}. "
            "Please run 'ppanggolin cluster' command first."
        )

    # Check if gene family sequences are available
    gene_family_sequences = pangenome.status.get("geneFamilySequences", "No")
    if gene_family_sequences == "No":
        raise AttributeError(
            f"No sequences associated with gene families for pangenome {getattr(pangenome, 'name', 'unknown')}. "
            "Please ensure gene family sequences are available."
        )


def write_clustering(clust_res: Path, outfile: Path) -> None:
    """
    Process and write clustering results with proper cluster IDs.

    This function reads raw MMseqs2 clustering results, assigns unique cluster IDs,
    and writes the processed results in a standardized format with proper headers.

    Args:
        clust_res: Path to the raw clustering results file from MMseqs2.
            Must be a TSV file with referent and member columns.
        outfile: Path where the processed clustering results will be written.
            Directory must exist and be writable.

    Raises:
        ClusteringError: If processing or writing fails.
        FileNotFoundError: If the input file doesn't exist.
    """
    # Validate input file
    if not clust_res.exists():
        raise FileNotFoundError(f"Clustering results file not found: {clust_res}")

    logger.debug(f"Processing clustering results from: {clust_res}")

    try:
        # Read raw clustering results
        # Expected format: referent \t member
        clust_df = pd.read_csv(
            clust_res,
            sep="\t",
            names=CONFIG.CLUSTER_COLUMN_NAMES[
                1:
            ],  # Skip the cluster_id column initially
            dtype=str,
        )

        if clust_df.empty:
            raise ClusteringError("Clustering results file is empty")

        logger.debug(f"Read {len(clust_df)} clustering entries")

        # Create cluster ID mapping
        # Each unique referent gets a unique cluster ID
        unique_referents = clust_df[CONFIG.CLUSTER_COLUMN_NAMES[1]].unique()
        cluster_id_mapping = {
            referent: cluster_id for cluster_id, referent in enumerate(unique_referents)
        }

        logger.debug(
            f"Created {len(cluster_id_mapping)} unique clusters",
        )

        # Create cluster ID dataframe
        cluster_id_df = pd.DataFrame(
            [
                {"cluster_id": cluster_id, "referent": referent}
                for referent, cluster_id in cluster_id_mapping.items()
            ]
        )

        # Merge cluster IDs with clustering results
        merged_clustering = cluster_id_df.merge(
            clust_df,
            on=CONFIG.CLUSTER_COLUMN_NAMES[1],  # referent column
            how="inner",
            validate="one_to_many",
        )

        # Ensure proper column order
        merged_clustering = merged_clustering[CONFIG.CLUSTER_COLUMN_NAMES]

        # Ensure output directory exists
        outfile.parent.mkdir(parents=True, exist_ok=True)

        # Write processed results
        merged_clustering.to_csv(outfile, sep="\t", header=True, index=False)

        total_entries = len(merged_clustering)
        unique_clusters = merged_clustering["cluster_id"].nunique()

        logger.info(
            "Clustering results processed: %d entries in %d clusters saved to: %s",
            total_entries,
            unique_clusters,
            outfile.absolute(),
        )

    except pd.errors.EmptyDataError:
        raise ClusteringError("Clustering results file is empty or malformed")
    except pd.errors.ParserError as e:
        raise ClusteringError(f"Error parsing clustering results file: {e}")
    except OSError as e:
        raise ClusteringError(f"File I/O error during clustering processing: {e}")
    except Exception as e:
        raise ClusteringError(f"Unexpected error during clustering processing: {e}")


def create_tsv(
    db: Path, clust: Path, output: Path, threads: int = CONFIG.DEFAULT_THREADS
) -> None:
    """
    Convert MMseqs2 clustering database to TSV format.

    This function uses MMseqs2's createtsv command to convert binary clustering
    results into a human-readable tab-separated values file.

    Args:
        db: Path to the MMseqs2 sequence database used for clustering.
        clust: Path to the MMseqs2 clustering results database.
        output: Path where the TSV file will be written.
        threads: Number of threads for conversion. Defaults to 1.

    Raises:
        ClusteringError: If TSV creation fails.
        FileNotFoundError: If input databases don't exist.
    """
    logger.debug("Converting clustering database to TSV format")

    try:
        # Prepare MMseqs2 createtsv command
        cmd = [
            "mmseqs",
            "createtsv",
            db.absolute().as_posix(),
            db.absolute().as_posix(),  # Use the same database for both query and target
            clust.absolute().as_posix(),
            output.absolute().as_posix(),
            "--threads",
            str(threads),
            "--full-header",  # Include full sequence headers
        ]

        logger.debug(f"Executing MMseqs2 createtsv: {' '.join(cmd)}")

        # Execute command with error handling
        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=False,
            text=True,
        )

        if result.returncode != 0:
            error_msg = f"MMseqs2 createtsv failed with return code {result.returncode}"
            if result.stderr:
                error_msg += f": {result.stderr.strip()}"
            raise ClusteringError(error_msg)

        logger.debug("TSV conversion completed successfully")

    except subprocess.SubprocessError as e:
        raise ClusteringError(f"Failed to execute MMseqs2 createtsv: {e}")
    except Exception as e:
        raise ClusteringError(f"Unexpected error during TSV creation: {e}")


def _execute_clustering_command(
    cmd: List[str],
    method_name: str,
) -> float:
    """
    Execute MMseqs2 clustering command with timing and error handling.

    Args:
        cmd: Command to execute as a list of strings.
        method_name: Name of the clustering method for logging.

    Returns:
        float: Execution time in seconds.

    Raises:
        ClusteringError: If command execution fails.
    """
    logger.debug(f"Executing {method_name} clustering: {' '.join(cmd)}")

    begin_time = time()
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=False,
            text=True,
        )

        if result.returncode != 0:
            error_msg = (
                f"MMseqs2 {method_name} failed with return code {result.returncode}"
            )
            if result.stderr:
                error_msg += f": {result.stderr.strip()}"
            raise ClusteringError(error_msg)

        execution_time = time() - begin_time

        return execution_time

    except subprocess.SubprocessError as e:
        raise ClusteringError(f"Failed to execute MMseqs2 {method_name}: {e}")


def linclust_launcher(
    seq_db: Path,
    mmseqs2_opt: Dict[str, Union[int, float, str]],
    lclust_db: Optional[Path] = None,
    tmpdir: Optional[Path] = None,
    threads: int = CONFIG.DEFAULT_THREADS,
) -> Path:
    """
    Launch MMseqs2 linclust (fast clustering) on gene family sequences.

    Linclust is a faster clustering method suitable for large datasets where
    speed is more important than sensitivity. It uses a linear time complexity
    algorithm for clustering.

    Args:
        seq_db:
            Path to MMseqs2 sequence database containing gene families.
        mmseqs2_opt:
            Dictionary containing MMseqs2 clustering parameters.
            Required keys: comp_bias_corr, kmer_per_seq, identity, coverage,
            cov_mode, eval, align_mode, max_seq_len, max_reject, clust_mode.
        lclust_db:
            Optional path for the clustering results database.
            If None, a temporary file will be created.
        tmpdir:
            Temporary directory for MMseqs2 operations.
            If None, a system temp directory will be used.
        threads:
            Number of threads for clustering. Defaults to 1.

    Returns:
        Path:
            Path to the MMseqs2 clustering results database.

    Raises:
        ClusteringError:
            If clustering execution fails.
        ClusteringValidationError:
            If parameters are invalid.
        FileNotFoundError:
            If the sequence database doesn't exist.
    """
    # Validate inputs
    if not seq_db.exists():
        raise FileNotFoundError(f"Sequence database not found: {seq_db}")

    _validate_clustering_parameters(threads, mmseqs2_opt)

    # Set the default temporary directory
    if tmpdir is None:
        tmpdir = Path(tempfile.gettempdir())
    else:
        _validate_directory_access(tmpdir, create_if_missing=True)

    logger.debug("Starting linclust clustering of gene families")

    try:
        # Create the clustering database if not provided
        temp_lclust_db = tmpdir / "linclust_db" if lclust_db is None else lclust_db
        return _execute_linclust(seq_db, mmseqs2_opt, temp_lclust_db, tmpdir, threads)

    except Exception as e:
        raise ClusteringError(f"Linclust clustering failed: {e}") from e


def _execute_linclust(
    seq_db: Path,
    mmseqs2_opt: Dict[str, Union[int, float, str]],
    lclust_db: Path,
    tmpdir: Path,
    threads: int,
) -> Path:
    """
    Execute the actual linclust command.

    Args:
        seq_db: Sequence database path.
        mmseqs2_opt: MMseqs2 options.
        lclust_db: Linclust database path.
        tmpdir: Temporary directory.
        threads: Number of threads.

    Returns:
        Path: Path to the clustering database.

    Raises:
        ClusteringError: If clustering execution fails.
    """
    # Prepare linclust command
    cmd = [
        "mmseqs",
        "cluster",
        str(seq_db.absolute()),
        str(lclust_db.absolute()),
        str(tmpdir.absolute()),
        "--threads",
        str(threads),
        "--comp-bias-corr",
        str(mmseqs2_opt["comp_bias_corr"]),
        "--kmer-per-seq",
        str(mmseqs2_opt["kmer_per_seq"]),
        "--min-seq-id",
        str(mmseqs2_opt["identity"]),
        "-c",
        str(mmseqs2_opt["coverage"]),
        "--cov-mode",
        str(mmseqs2_opt["cov_mode"]),
        "-e",
        str(mmseqs2_opt["eval"]),
        "--alignment-mode",
        str(mmseqs2_opt["align_mode"]),
        "--max-seq-len",
        str(mmseqs2_opt["max_seq_len"]),
        "--max-rejected",
        str(mmseqs2_opt["max_reject"]),
        "--cluster-mode",
        str(mmseqs2_opt["clust_mode"]),
    ]

    # Execute clustering
    execution_time = _execute_clustering_command(cmd, "linclust")

    logger.info("Linclust clustering completed in %.2f seconds", execution_time)
    return lclust_db


def cluster_launcher(
    seq_db: Path,
    mmseqs2_opt: Dict[str, Union[int, float, str]],
    cluster_db: Optional[Path] = None,
    tmpdir: Optional[Path] = None,
    threads: int = CONFIG.DEFAULT_THREADS,
) -> Path:
    """
    Launch MMseqs2 cluster (sensitive clustering) on gene family sequences.

    Cluster is a more sensitive but slower clustering method that provides
    better clustering quality through more thorough sequence comparison.

    Args:
        seq_db:
            Path to MMseqs2 sequence database containing gene families.
        mmseqs2_opt:
            Dictionary containing MMseqs2 clustering parameters. Required keys:
            - max_seqs,
            - min_ungapped,
            - comp_bias_corr,
            - sensitivity,
            - kmer_per_seq,
            - identity,
            - coverage,
            - cov_mode,
            - eval,
            - align_mode,
            - max_seq_len,
            - max_reject,
            - clust_mode.
        cluster_db:
            Optional path for the clustering results database.
            If None, a temporary file will be created.
        tmpdir:
            Temporary directory for MMseqs2 operations.
            If None, the system temp directory will be used.
        threads:
            Number of threads for clustering. Defaults to 1.

    Returns:
        Path:
            Path to the MMseqs2 clustering results database.

    Raises:
        ClusteringError:
            If clustering execution fails.
        ClusteringValidationError:
            If parameters are invalid.
        FileNotFoundError:
            If the sequence database doesn't exist.
    """
    # Validate inputs
    if not seq_db.exists():
        raise FileNotFoundError(f"Sequence database not found: {seq_db}")

    _validate_clustering_parameters(threads, mmseqs2_opt)

    # Validate cluster-specific parameters
    cluster_specific = ["max_seqs", "min_ungapped", "sensitivity"]
    for param in cluster_specific:
        if param not in mmseqs2_opt:
            raise ClusteringValidationError(
                f"Missing cluster-specific parameter: {param}"
            )

    # Set the default temporary directory
    if tmpdir is None:
        tmpdir = Path(tempfile.gettempdir())
    else:
        _validate_directory_access(tmpdir, create_if_missing=True)

    logger.debug("Starting sensitive cluster clustering of gene families")

    try:
        # Create the clustering database if not provided
        temp_cluster_db = tmpdir / "cluster_db" if cluster_db is None else cluster_db
        return _execute_cluster(seq_db, mmseqs2_opt, temp_cluster_db, tmpdir, threads)

    except Exception as e:
        raise ClusteringError(f"Cluster clustering failed: {e}") from e


def _execute_cluster(
    seq_db: Path,
    mmseqs2_opt: Dict[str, Union[int, float, str]],
    cluster_db: Path,
    tmpdir: Path,
    threads: int,
) -> Path:
    """
    Execute the actual cluster command.

    Args:
        seq_db: Sequence database path.
        mmseqs2_opt: MMseqs2 options.
        cluster_db: Cluster database path.
        tmpdir: Temporary directory.
        threads: Number of threads.

    Returns:
        Path: Path to the clustering database.

    Raises:
        ClusteringError: If clustering execution fails.
    """
    # Prepare cluster command
    cmd = [
        "mmseqs",
        "cluster",
        str(seq_db.absolute()),
        str(cluster_db.absolute()),
        str(tmpdir.absolute()),
        "--threads",
        str(threads),
        "--max-seqs",
        str(mmseqs2_opt["max_seqs"]),
        "--min-ungapped-score",
        str(mmseqs2_opt["min_ungapped"]),
        "--comp-bias-corr",
        str(mmseqs2_opt["comp_bias_corr"]),
        "-s",
        str(mmseqs2_opt["sensitivity"]),
        "--kmer-per-seq",
        str(mmseqs2_opt["kmer_per_seq"]),
        "--min-seq-id",
        str(mmseqs2_opt["identity"]),
        "-c",
        str(mmseqs2_opt["coverage"]),
        "--cov-mode",
        str(mmseqs2_opt["cov_mode"]),
        "-e",
        str(mmseqs2_opt["eval"]),
        "--alignment-mode",
        str(mmseqs2_opt["align_mode"]),
        "--max-seq-len",
        str(mmseqs2_opt["max_seq_len"]),
        "--max-rejected",
        str(mmseqs2_opt["max_reject"]),
        "--cluster-mode",
        str(mmseqs2_opt["clust_mode"]),
    ]

    # Execute clustering
    execution_time = _execute_clustering_command(cmd, "cluster")

    logger.info("Cluster clustering completed in %.2f seconds", execution_time)
    return cluster_db


def cluster_gene_families(
    pangenomes: Pangenomes,
    method: str,
    mmseqs2_opt: Dict[str, Union[int, float, str]],
    tmpdir: Path,
    threads: int = CONFIG.DEFAULT_THREADS,
    lock: Optional[Lock] = None,
    disable_bar: bool = False,
) -> Path:
    """
    Cluster gene families from multiple pangenomes using MMseqs2.

    This is the main function that orchestrates the complete clustering workflow:
    writing sequences, creating databases, performing clustering, and formatting results.

    Args:
        pangenomes:
            Pangenomes object containing multiple pangenome instances.
            All pangenomes must have clustered genes and family sequences.
        method:
            Clustering method to use. Must be "linclust" or "cluster".
            - "linclust" is faster but less sensitive,
            - "cluster" is more sensitive but slower.
        mmseqs2_opt:
            Dictionary containing MMseqs2 clustering parameters.
            Must include all required parameters for the chosen method.
        tmpdir:
            Temporary directory for operations. If None, uses system temp.
        threads:
            Number of threads for processing. Defaults to 1.
        lock:
            Optional multiprocessing Lock for thread safety. If None, operations
            may not be thread-safe in multiprocessing contexts.
        disable_bar:
            Whether to disable progress bars. Defaults to False.

    Returns:
        Path:
            Path to the final clustering results file in TSV format with columns:
            cluster_id, referent, in_clust.

    Raises:
        ClusteringError:
            If the clustering process fails.
        ClusteringValidationError:
            If parameters are invalid.
        ValueError:
            If the method is not "linclust" or "cluster".

    Notes:
        Temporary directory is supposed to be already validated
    """
    if method not in ClusteringMethod.CHOICES:
        raise ValueError(
            f"Invalid clustering method: {method}. "
            f"Must choose between '{ClusteringMethod.LINCLUST}' or '{ClusteringMethod.CLUSTER}'. "
            "See MMseqs2 documentation for more information."
        )

    _validate_clustering_parameters(threads, mmseqs2_opt)

    logger.info(f"Starting gene family clustering with method: {method}")
    logger.info(f"Processing {len(pangenomes)} pangenomes")

    try:
        # Write pangenome families sequences
        logger.info("Writing pangenome families sequences...")
        pangenome2families_seq = write_pangenomes_families_sequences(
            pangenomes=pangenomes,
            output=tmpdir,
            threads=threads,
            lock=lock,
            disable_bar=disable_bar,
        )

        logger.info(
            f"Successfully wrote sequences for {len(pangenome2families_seq)} pangenomes"
        )

        # Create combined MMseqs2 database
        logger.info("Creating combined MMseqs2 database...")
        sequence_files = list(pangenome2families_seq.values())

        try:
            merged_db = createdb(seq_files=sequence_files, output=tmpdir)
            logger.debug(f"Combined database created: {merged_db}")

        except (PangenomeProcessingError, MMSeqsError) as e:
            raise ClusteringError(f"Failed to create sequence database: {e}")

        # Perform clustering based on the selected method
        logger.info("Beginning gene family clustering...")

        if method == ClusteringMethod.LINCLUST:
            cluster_db = linclust_launcher(
                seq_db=merged_db,
                mmseqs2_opt=mmseqs2_opt,
                tmpdir=tmpdir,
                threads=threads,
            )
        elif method == ClusteringMethod.CLUSTER:
            cluster_db = cluster_launcher(
                seq_db=merged_db,
                mmseqs2_opt=mmseqs2_opt,
                tmpdir=tmpdir,
                threads=threads,
            )
        else:
            # This should not happen due to earlier validation, but included for completeness
            raise ClusteringError(f"Unsupported clustering method: {method}")

        # Convert clustering results to TSV format
        logger.debug("Converting clustering results to TSV format...")
        cluster_results_file = tmpdir / "cluster_results.tsv"
        create_tsv(
            db=merged_db, clust=cluster_db, output=cluster_results_file, threads=threads
        )

        logger.info("Gene family clustering completed successfully")
        return cluster_results_file

    except Exception as e:
        raise ClusteringError(f"Gene family clustering process failed: {e}") from e


def launch(args: argparse.Namespace) -> None:
    """
    Main entry point for clustering operations.

    This function orchestrates the complete clustering the workflow, including
    parameter validation, pangenome loading, clustering execution, and result
    processing based on the specified method and parameters.

    Args:
        args: Parsed command line arguments containing all configuration
            parameters for the clustering operation.
            Expected attributes: pangenomes, output, method, tmpdir, threads,
            keep_tmp, disable_prog_bar, force, and MMseqs2 parameters.

    Raises:
        ClusteringError: If any step of the clustering process fails.
        ClusteringValidationError: If parameter validation fails.
        NotADirectoryError: If specified directories are invalid.
    """
    # Validate command line parameters
    logger.info("Validating clustering parameters...")
    check_cluster_parameters(args)

    # Create output directory
    logger.debug(f"Creating output directory: {args.output}")
    mkdir(args.output, getattr(args, "force", False))

    # Set up multiprocessing
    manager = Manager()
    lock = manager.Lock()

    # Define requirements for pangenome loading
    need_info = {"need_families": True, "need_families_sequences": False}

    # Load pangenomes with validation
    logger.info(f"Loading pangenomes from: {args.pangenomes}")
    pangenomes = load_pangenomes(
        pangenome_list=args.pangenomes,
        need_info=need_info,
        check_function=check_pangenome_cluster,
        max_workers=getattr(args, "threads", CONFIG.DEFAULT_THREADS),
        lock=lock,
        disable_bar=getattr(args, "disable_prog_bar", False),
    )

    logger.info(f"Successfully loaded {len(pangenomes)} pangenomes")

    # Create the temporary directory
    tmpdir = Path(tempfile.mkdtemp(dir=args.tmpdir))
    logger.debug(f"Created temporary directory: {tmpdir}")

    try:
        # Prepare MMseqs2 options from command line arguments
        mmseqs2_options = _prepare_mmseqs2_options(args)

        # Perform clustering
        logger.info(f"Starting clustering with method: {args.method}")
        clustering_results = cluster_gene_families(
            pangenomes=pangenomes,
            method=args.method,
            mmseqs2_opt=mmseqs2_options,
            tmpdir=tmpdir,
            threads=getattr(args, "threads", CONFIG.DEFAULT_THREADS),
            lock=lock,
            disable_bar=getattr(args, "disable_prog_bar", False),
        )

        # Process and write final results
        logger.info("Processing final clustering results...")
        final_output = args.output / CONFIG.CLUSTERING_OUTPUT
        write_clustering(clustering_results, final_output)

        logger.info("Clustering process completed successfully")
    except ClusteringError as e:
        raise ClusteringError(f"Clustering process failed: {e}") from e
    except Exception as e:
        raise ClusteringError(
            f"Clustering process failed with unexpected error: {e}"
        ) from e
    finally:
        # Clean up the temporary directory unless requested to keep it
        if not args.keep_tmp:
            logger.debug(f"Cleaning up temporary directory: {tmpdir}")
            rmtree(tmpdir, ignore_errors=True)
        else:
            logger.info(f"Temporary files kept in: {tmpdir}")


def _prepare_mmseqs2_options(
    args: argparse.Namespace,
) -> Dict[str, Union[int, float, str]]:
    """
    Prepare MMseqs2 options dictionary from command line arguments.

    Args:
        args: Command line arguments containing MMseqs2 parameters.

    Returns:
        Dict[str, Union[int, float, str]]: Dictionary of MMseqs2 options.
    """
    # Start with default options
    mmseqs2_options = CONFIG.DEFAULT_MMSEQS2_OPTIONS.copy()

    # Override with command line arguments if provided
    option_mappings = {
        "cluster_identity": "identity",
        "cluster_coverage": "coverage",
        "cluster_cov_mode": "cov_mode",
        "cluster_eval": "eval",
        "cluster_sensitivity": "sensitivity",
        "cluster_max_seqs": "max_seqs",
        "cluster_comp_bias_corr": "comp_bias_corr",
        "cluster_kmer_per_seq": "kmer_per_seq",
        "cluster_align_mode": "align_mode",
        "cluster_max_seq_len": "max_seq_len",
        "cluster_max_reject": "max_reject",
        "cluster_clust_mode": "clust_mode",
        "cluster_min_ungapped": "min_ungapped",
    }

    for arg_name, option_key in option_mappings.items():
        if hasattr(args, arg_name):
            value = getattr(args, arg_name)
            if value is not None:
                mmseqs2_options[option_key] = value

    return mmseqs2_options


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Create the argument subparser for clustering command.

    This function sets up the command-line interface for the clustering
    functionality within the PANORAMA tool suite.

    Args:
        sub_parser: The subparser object from argparse to add clustering command to.

    Returns:
        argparse.ArgumentParser: Configured parser for clustering command.
    """
    parser = sub_parser.add_parser(
        "cluster",
        help="Cluster gene families across pangenomes using MMseqs2",
        description="Perform gene family clustering across multiple pangenomes "
        "using MMseqs2 with support for both fast (linclust) and "
        "sensitive (cluster) clustering methods.",
    )
    parser_cluster(parser)
    return parser


def parser_mmseqs2_cluster(parser: argparse.ArgumentParser) -> argparse._ArgumentGroup:
    """
    Add MMseqs2-specific clustering arguments to the parser.

    Args:
        parser: ArgumentParser to add MMseqs2 arguments to.

    Returns:
        argparse._ArgumentGroup: The argument group containing MMseqs2 clustering options.
    """
    mmseqs_group = parser.add_argument_group(
        title="MMseqs2 clustering parameters",
        description="Configure MMseqs2 clustering behavior. "
        "See MMseqs2 documentation for detailed parameter descriptions.",
    )

    # Core clustering parameters
    mmseqs_group.add_argument(
        "--cluster_identity",
        required=False,
        type=float,
        default=CONFIG.DEFAULT_IDENTITY,
        metavar="FLOAT",
        help=f"Minimum sequence identity threshold (0.0-1.0). Default: {CONFIG.DEFAULT_IDENTITY}",
    )

    mmseqs_group.add_argument(
        "--cluster_coverage",
        required=False,
        type=float,
        default=CONFIG.DEFAULT_COVERAGE,
        metavar="FLOAT",
        help=f"Minimum coverage threshold (0.0-1.0). Default: {CONFIG.DEFAULT_COVERAGE}",
    )

    mmseqs_group.add_argument(
        "--cluster_cov_mode",
        required=False,
        type=int,
        default=CONFIG.DEFAULT_COV_MODE,
        choices=[0, 1, 2, 3, 4, 5],
        metavar="INT",
        help=f"Coverage mode: 0=query, 1=target, 2=shorter seq, 3=longer seq, "
        f"4=query and target, 5=shorter and longer seq. Default: {CONFIG.DEFAULT_COV_MODE}",
    )

    mmseqs_group.add_argument(
        "--cluster_eval",
        required=False,
        type=float,
        default=CONFIG.DEFAULT_EVAL,
        metavar="FLOAT",
        help=f"E-value threshold. Default: {CONFIG.DEFAULT_EVAL}",
    )

    # Sensitivity and performance parameters
    mmseqs_group.add_argument(
        "--cluster_sensitivity",
        required=False,
        type=float,
        metavar="FLOAT",
        help="Search sensitivity (cluster method only). Higher values = more sensitive but slower",
    )

    mmseqs_group.add_argument(
        "--cluster_max_seqs",
        required=False,
        type=int,
        metavar="INT",
        help="Maximum number of sequences per cluster representative (cluster method only)",
    )

    # Advanced parameters
    mmseqs_group.add_argument(
        "--cluster_comp_bias_corr",
        required=False,
        type=int,
        choices=[0, 1],
        metavar="INT",
        help="Compositional bias correction: 0=disabled, 1=enabled",
    )

    mmseqs_group.add_argument(
        "--cluster_kmer_per_seq",
        required=False,
        type=int,
        metavar="INT",
        help="Number of k-mers per sequence",
    )

    mmseqs_group.add_argument(
        "--cluster_align_mode",
        required=False,
        type=int,
        choices=[0, 1, 2, 3, 4],
        metavar="INT",
        help="Alignment mode: 0=automatic, 1=only score, 2=only extended, 3=score+extended, 4=fast+extended",
    )

    mmseqs_group.add_argument(
        "--cluster_max_seq_len",
        required=False,
        type=int,
        metavar="INT",
        help="Maximum sequence length",
    )

    mmseqs_group.add_argument(
        "--cluster_max_reject",
        required=False,
        type=int,
        metavar="INT",
        help="Maximum number of rejected sequences",
    )

    mmseqs_group.add_argument(
        "--cluster_clust_mode",
        required=False,
        type=int,
        choices=[0, 1, 2, 3],
        metavar="INT",
        help="Clustering mode: 0=Set Cover, 1=Connected Component, 2=Greedy, 3=Greedy Low Memory",
    )

    mmseqs_group.add_argument(
        "--cluster_min_ungapped",
        required=False,
        type=int,
        metavar="INT",
        help="Minimum ungapped alignment score (cluster method only)",
    )

    return mmseqs_group


def parser_cluster(parser: argparse.ArgumentParser) -> None:
    """
    Configure the argument parser for clustering command.

    This function adds all necessary command-line arguments for the clustering
    functionality, including required arguments, clustering methods, and optional
    parameters.

    Args:
        parser: ArgumentParser to configure with clustering arguments.
    """
    # Required arguments group
    required_group = parser.add_argument_group(
        title="Required arguments",
        description="All of the following arguments are required:",
    )

    required_group.add_argument(
        "-p",
        "--pangenomes",
        required=True,
        type=Path,
        metavar="FILE",
        help="Path to TSV file containing list of pangenome .h5 files to process",
    )

    required_group.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        metavar="DIR",
        help="Output directory where clustering results will be written",
    )

    required_group.add_argument(
        "-m",
        "--method",
        required=True,
        choices=ClusteringMethod.CHOICES,
        metavar="METHOD",
        help=f"Clustering method: {ClusteringMethod.LINCLUST} (fast) or "
        f"{ClusteringMethod.CLUSTER} (sensitive)",
    )

    # Add MMseqs2 specific arguments
    parser_mmseqs2_cluster(parser)

    # Optional arguments group
    optional_group = parser.add_argument_group(title="Optional arguments")

    optional_group.add_argument(
        "--tmpdir",
        required=False,
        type=Path,
        default=Path(tempfile.gettempdir()),
        metavar="DIR",
        help=f"Directory for temporary files. Default: {tempfile.gettempdir()}",
    )

    optional_group.add_argument(
        "--keep_tmp",
        required=False,
        default=False,
        action="store_true",
        help="Keep temporary files after completion (useful for debugging)",
    )

    optional_group.add_argument(
        "--threads",
        required=False,
        type=int,
        default=CONFIG.DEFAULT_THREADS,
        metavar="INT",
        help=f"Number of threads for parallel processing. Default: {CONFIG.DEFAULT_THREADS}",
    )
