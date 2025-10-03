#!/usr/bin/env python3
# coding:utf-8

"""
PANORAMA alignment module for inter-pangenome gene family comparisons.

This module provides functionality to align gene families between pangenomes using
MMseqs2, supporting both inter-pangenome-only and all-against-all alignment modes.
It handles sequence database creation, alignment execution, and result processing
with comprehensive error handling and progress tracking.
"""

# default libraries
from __future__ import annotations
import argparse
import logging
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from itertools import combinations
from shutil import rmtree
from multiprocessing import Manager, Lock
from time import time
from dataclasses import dataclass

# installed libraries
from tqdm import tqdm
import pandas as pd

# local libraries
from panorama.utils import mkdir, init_lock
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.format.read_binaries import load_pangenomes
from panorama.alignment.utils import (
    createdb,
    write_pangenomes_families_sequences,
    PangenomeProcessingError,
    MMSeqsError,
)

logger = logging.getLogger("PANORAMA")


# Configuration constants
@dataclass(frozen=True)
class AlignmentConfig:
    """Configuration constants for alignment operations."""

    # MMseqs2 output format for alignment results
    ALIGN_FORMAT: Optional[List[str]] = None
    ALIGN_COLUMNS: Optional[List[str]] = None

    # Default alignment parameters
    DEFAULT_IDENTITY: float = 0.8
    DEFAULT_COVERAGE: float = 0.8
    DEFAULT_COV_MODE: int = 0
    DEFAULT_THREADS: int = 1

    # Output filenames
    INTER_PANGENOMES_OUTPUT: str = "inter_pangenomes.tsv"
    ALL_AGAINST_ALL_OUTPUT: str = "all_against_all.tsv"

    def __post_init__(self):
        # Initialize mutable defaults
        if self.ALIGN_FORMAT is None:
            object.__setattr__(
                self,
                "ALIGN_FORMAT",
                [
                    "query",
                    "target",
                    "fident",
                    "qlen",
                    "tlen",
                    "alnlen",
                    "evalue",
                    "bits",
                ],
            )
        if self.ALIGN_COLUMNS is None:
            object.__setattr__(
                self,
                "ALIGN_COLUMNS",
                [
                    "query",
                    "target",
                    "identity",
                    "qlength",
                    "tlength",
                    "alnlength",
                    "e_value",
                    "bits",
                ],
            )


# Global configuration instance
CONFIG = AlignmentConfig()


class AlignmentError(Exception):
    """Custom exception for alignment-related errors."""

    pass


class AlignmentValidationError(AlignmentError):
    """Custom exception for alignment parameter validation errors."""

    pass


def _validate_alignment_parameters(
    identity: float, coverage: float, cov_mode: int, threads: int
) -> None:
    """
    Validate alignment parameters.

    Args:
        identity: Sequence identity threshold (0.0-1.0).
        coverage: Coverage threshold (0.0-1.0).
        cov_mode: Coverage mode for MMseqs2 (0-5).
        threads: Number of threads (positive integer).

    Raises:
        AlignmentValidationError: If any parameter is invalid.
    """
    if not isinstance(identity, (int, float)) or not (0.0 <= identity <= 1.0):
        raise AlignmentValidationError(
            f"Identity must be between 0.0 and 1.0, got: {identity}"
        )

    if not isinstance(coverage, (int, float)) or not (0.0 <= coverage <= 1.0):
        raise AlignmentValidationError(
            f"Coverage must be between 0.0 and 1.0, got: {coverage}"
        )

    if not isinstance(cov_mode, int) or not (0 <= cov_mode <= 5):
        raise AlignmentValidationError(
            f"Coverage mode must be between 0 and 5, got: {cov_mode}"
        )

    if not isinstance(threads, int) or threads <= 0:
        raise AlignmentValidationError(
            f"Threads must be positive integer, got: {threads}"
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
        AlignmentValidationError: If directory validation fails.
    """
    if not isinstance(directory, Path):
        raise AlignmentValidationError(f"Expected Path object, got: {type(directory)}")

    if not directory.exists():
        if create_if_missing:
            try:
                directory.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                raise AlignmentValidationError(
                    f"Cannot create directory {directory}: {e}"
                )
        else:
            raise AlignmentValidationError(f"Directory does not exist: {directory}")

    if not directory.is_dir():
        raise AlignmentValidationError(f"Path is not a directory: {directory}")


def check_align_parameters(args: argparse.Namespace) -> None:
    """
    Validate command line arguments for alignment operations.

    This function performs comprehensive validation of all alignment parameters
    provided via command line arguments, ensuring they meet the required
    constraints before proceeding with alignment operations.

    Args:
        args: Parsed command line arguments containing alignment parameters.
            Expected attributes: tmpdir, align_identity, align_coverage.

    Raises:
        AlignmentValidationError: If any parameter validation fails.
        NotADirectoryError: If tmpdir is not a valid directory.
    """

    # Validate temporary directory
    if not hasattr(args, "tmpdir") or not isinstance(args.tmpdir, Path):
        raise AlignmentValidationError("tmpdir must be a valid Path object")

    try:
        _validate_directory_access(args.tmpdir, create_if_missing=True)
    except AlignmentValidationError as e:
        raise NotADirectoryError(
            f"The given path for temporary directory is not valid: {args.tmpdir}. "
            f"Error: {e}. Please check your path and try again."
        )

    # Validate alignment parameters
    identity = getattr(args, "align_identity", CONFIG.DEFAULT_IDENTITY)
    coverage = getattr(args, "align_coverage", CONFIG.DEFAULT_COVERAGE)

    if identity > 1 or coverage > 1:
        raise AlignmentValidationError(
            f"Identity ({identity}) and coverage ({coverage}) must be between 0 and 1"
        )

    logger.debug("Alignment parameters validated successfully")


def check_pangenome_align(pangenome: Pangenome) -> None:
    """
    Validate that a pangenome is ready for alignment operations.

    This function checks that the pangenome has been properly processed
    and contains the necessary data for gene family alignment operations.

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


def write_alignment(
    query_db: Path,
    target_db: Path,
    aln_db: Path,
    outfile: Path,
    threads: int = CONFIG.DEFAULT_THREADS,
) -> None:
    """
    Convert MMseqs2 alignment database to human-readable format.

    This function uses MMseqs2's convertalis command to convert binary
    alignment results into a tab-separated values file with specified
    output format for downstream analysis.

    Args:
        query_db: Path to MMseqs2 query database.
        target_db: Path to MMseqs2 target database.
        aln_db: Path to MMseqs2 alignment results database.
        outfile: Path where the converted alignment results will be written.
        threads: Number of threads for the conversion process. Defaults to 1.

    Raises:
        AlignmentError: If alignment conversion fails.
        FileNotFoundError: If input databases don't exist.
    """

    # Validate parameters
    try:
        # Prepare MMseqs2 convertalis command
        cmd = [
            "mmseqs",
            "convertalis",
            str(query_db.absolute()),
            str(target_db.absolute()),
            str(aln_db.absolute()),
            str(outfile.absolute()),
            "--format-output",
            ",".join(CONFIG.ALIGN_FORMAT),
            "--threads",
            str(threads),
        ]

        logger.debug(f"Executing alignment conversion: {' '.join(cmd)}")
        logger.info("Extracting alignments...")

        # Execute command with error handling
        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=False,
            text=True,
        )

        if result.returncode != 0:
            error_msg = (
                f"MMseqs2 convertalis failed with return code {result.returncode}"
            )
            if result.stderr:
                error_msg += f": {result.stderr.strip()}"
            raise AlignmentError(error_msg)

        logger.debug("Alignment conversion completed successfully")

    except subprocess.SubprocessError as e:
        raise AlignmentError(f"Failed to execute MMseqs2 convertalis: {e}")
    except Exception as e:
        raise AlignmentError(f"Unexpected error during alignment conversion: {e}")


def _execute_alignment(
    query_db: Path,
    target_db: Path,
    aln_db: Path,
    tmpdir: Path,
    identity: float,
    coverage: float,
    cov_mode: int,
    threads: int,
) -> Path:
    """
    Execute the actual MMseqs2 alignment command.

    Args:
        query_db: Query database path.
        target_db: Target database path.
        aln_db: Alignment database path.
        tmpdir: Temporary directory.
        identity: Identity threshold.
        coverage: Coverage threshold.
        cov_mode: Coverage mode.
        threads: Number of threads.

    Returns:
        Path: Path to the alignment database.

    Raises:
        AlignmentError: If alignment execution fails.
    """
    # Prepare MMseqs2 search command
    cmd = [
        "mmseqs",
        "search",
        query_db.absolute().as_posix(),
        target_db.absolute().as_posix(),
        aln_db.absolute().as_posix(),
        tmpdir.absolute().as_posix(),
        "-a",  # Include backtrace information
        "--min-seq-id",
        str(identity),
        "-c",
        str(coverage),
        "--cov-mode",
        str(cov_mode),
        "--threads",
        str(threads),
    ]
    logger.debug(f"Executing MMseqs2 search: {' '.join(cmd)}")

    # Execute alignment with timing
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
            error_msg = f"MMseqs2 search failed with return code {result.returncode}"
            if result.stderr:
                error_msg += f": {result.stderr.strip()}"
            raise AlignmentError(error_msg)

        align_time = time() - begin_time
        logger.debug("Alignment completed in %.2f seconds", align_time)

        return aln_db

    except subprocess.SubprocessError as e:
        raise AlignmentError(f"Failed to execute MMseqs2 search: {e}")


def align_db(
    query_db: Path,
    target_db: Path,
    tmpdir: Path,
    aln_db: Optional[Path] = None,
    query_name: str = "",
    target_name: str = "",
    identity: float = CONFIG.DEFAULT_IDENTITY,
    coverage: float = CONFIG.DEFAULT_COVERAGE,
    cov_mode: int = CONFIG.DEFAULT_COV_MODE,
    threads: int = CONFIG.DEFAULT_THREADS,
) -> Path:
    """
    Perform sequence alignment between query and target databases using MMseqs2.

    This function executes MMseqs2 search to align sequences from the query
    database against the target database, applying specified identity and
    coverage thresholds.

    Args:
        query_name:
            Name of the query database
        target_name:
            Name of the target database
        query_db:
            Path to MMseqs2 query sequences database.
        target_db:
            Path to MMseqs2 target sequences database.
        aln_db:
            Optional path for the alignment results database. If None, a temporary file will be created.
        tmpdir:
            Temporary directory for MMseqs2 operations.
            If None, the system temp directory will be used.
        identity:
            Minimum sequence identity threshold (0.0-1.0). Defaults to 0.8.
        coverage:
            Minimum coverage threshold (0.0-1.0). Defaults to 0.8.
        cov_mode:
            Coverage mode for MMseqs2 (0-5). Defaults to 0.
                - 0: coverage of the query,
                - 1: coverage of the target,
                - 2: coverage of the shorter sequence.
        threads:
            Number of threads for alignment. Defaults to 1.

    Returns:
        Path: Path to the MMseqs2 alignment results database.

    Raises:
        AlignmentError: If alignment execution fails.
        AlignmentValidationError: If parameters are invalid.
        FileNotFoundError: If input databases don't exist.
    TODO:
        - see how validate input databases
    """

    # Validate input parameters
    _validate_alignment_parameters(identity, coverage, cov_mode, threads)

    logger.debug(f"Starting database alignment between {query_db} and {target_db}")

    try:
        # Create the alignment database if not provided
        if aln_db is None:
            temp_aln_db_name = "align_" + (
                f"{query_name}_vs_{target_name}"
                if query_name and target_name
                else query_name or target_name or ""
            ).rstrip("_")
            temp_aln_db = tmpdir / temp_aln_db_name
            return _execute_alignment(
                query_db,
                target_db,
                temp_aln_db,
                tmpdir,
                identity,
                coverage,
                cov_mode,
                threads,
            )
        else:
            return _execute_alignment(
                query_db,
                target_db,
                aln_db,
                tmpdir,
                identity,
                coverage,
                cov_mode,
                threads,
            )

    except Exception as e:
        raise AlignmentError(f"Database alignment failed: {e}") from e


def align_pangenomes_pair(
    pangenomes_pair: Tuple[str, str],
    tmpdir: Path,
    db_pair: Tuple[Path, Path],
    identity: float = CONFIG.DEFAULT_IDENTITY,
    coverage: float = CONFIG.DEFAULT_COVERAGE,
    cov_mode: int = CONFIG.DEFAULT_COV_MODE,
    threads: int = CONFIG.DEFAULT_THREADS,
) -> Path:
    """
    Align gene families between two specific pangenomes.

    This function performs pairwise alignment between gene families from
    two different pangenomes, executing both the alignment and conversion
    to human-readable format.

    Args:
        pangenomes_pair: Tuple containing names of the two pangenomes to align.
        db_pair: Tuple containing paths to MMseqs2 databases for the pangenome pair.
        identity: Minimum sequence identity threshold (0.0-1.0). Defaults to 0.8.
        coverage: Minimum coverage threshold (0.0-1.0). Defaults to 0.8.
        cov_mode: Coverage mode for MMseqs2 (0-5). Defaults to 0.
        tmpdir: Temporary directory for operations. If None, uses system temp.
        threads: Number of threads for processing. Defaults to 1.

    Returns:
        Path: Path to the alignment results file in TSV format.

    Raises:
        AlignmentError: If alignment or conversion fails.
        AlignmentValidationError: If parameters are invalid.
    """

    # Validate input parameters
    _validate_alignment_parameters(identity, coverage, cov_mode, threads)

    if len(pangenomes_pair) != 2 or len(db_pair) != 2:
        raise AlignmentValidationError(
            "pangenomes_pair and db_pair must contain exactly 2 elements"
        )

    pangenome1, pangenome2 = pangenomes_pair
    query_db, target_db = db_pair

    logger.debug(f"Aligning gene families between {pangenome1} and {pangenome2}")
    try:
        # Perform alignment
        aln_db = align_db(
            query_db=query_db,
            target_db=target_db,
            tmpdir=tmpdir,
            query_name=pangenome1,
            target_name=pangenome2,
            identity=identity,
            coverage=coverage,
            cov_mode=cov_mode,
            threads=threads,
        )

        logger.debug(f"Writing alignment results between {pangenome1} and {pangenome2}")

        # Create the output file for alignment results
        aln_results_file = tmpdir / f"{pangenome1}_vs_{pangenome2}.tsv"
        write_alignment(query_db, target_db, aln_db, aln_results_file, threads)

        logger.debug(f"Alignment between {pangenome1} and {pangenome2} completed")
        return aln_results_file

    except Exception as e:
        raise AlignmentError(
            f"Pairwise alignment failed for {pangenome1} vs {pangenome2}: {e}"
        ) from e


def align_pangenomes(
    pangenome2db: Dict[str, Path],
    tmpdir: Path,
    identity: float = CONFIG.DEFAULT_IDENTITY,
    coverage: float = CONFIG.DEFAULT_COVERAGE,
    cov_mode: int = CONFIG.DEFAULT_COV_MODE,
    threads: int = CONFIG.DEFAULT_THREADS,
    disable_bar: bool = False,
) -> List[Path]:
    """
    Perform all pairwise alignments between multiple pangenomes.

    This function executes alignments between all possible pairs of pangenomes
    using their MMseqs2 databases, with progress tracking and error handling.

    Args:
        pangenome2db: Dictionary mapping pangenome names to their MMseqs2 database paths.
        identity: Minimum sequence identity threshold (0.0-1.0). Defaults to 0.8.
        coverage: Minimum coverage threshold (0.0-1.0). Defaults to 0.8.
        cov_mode: Coverage mode for MMseqs2 (0-5). Defaults to 0.
        tmpdir: Temporary directory for operations. If None, uses system temp.
        threads: Number of threads per alignment. Defaults to 1.
        disable_bar: Whether to disable the progress bar. Defaults to False.

    Returns:
        List[Path]: List of paths to alignment result files for each pangenome pair.

    Raises:
        AlignmentError: If any alignment fails.
        AlignmentValidationError: If parameters are invalid.
    """
    # Validate input parameters
    if len(pangenome2db) < 2:
        raise AlignmentValidationError(
            "At least 2 pangenomes required for pairwise alignment"
        )

    _validate_alignment_parameters(identity, coverage, cov_mode, threads)

    # Generate all possible pangenome pairs
    pangenomes_pairs = list(combinations(pangenome2db.keys(), 2))
    total_pairs = len(pangenomes_pairs)

    logger.info(
        f"Aligning gene families between {total_pairs} pangenome pairs with {threads} threads..."
    )

    results = []

    for pangenomes_pair in tqdm(
        pangenomes_pairs,
        unit="pangenome pairs",
        disable=disable_bar,
        desc="Aligning pangenomes",
    ):
        try:
            # Get database paths for this pair
            db_pair = (
                pangenome2db[pangenomes_pair[0]],
                pangenome2db[pangenomes_pair[1]],
            )

            # Perform pairwise alignment
            alignment_result = align_pangenomes_pair(
                pangenomes_pair=pangenomes_pair,
                tmpdir=tmpdir,
                db_pair=db_pair,
                identity=identity,
                coverage=coverage,
                cov_mode=cov_mode,
                threads=threads,
            )

            results.append(alignment_result)
            logger.debug(
                f"Completed alignment for pair: {pangenomes_pair[0]} vs {pangenomes_pair[1]}"
            )

        except Exception as e:
            error_msg = (
                f"Failed alignment for {pangenomes_pair[0]} vs {pangenomes_pair[1]}"
            )
            raise AlignmentError(error_msg) from e

    return results


def merge_aln_res(align_results: List[Path], outfile: Path) -> None:
    """
    Merge multiple alignment result files into a single consolidated file.

    This function reads multiple TSV alignment files and combines them into
    a single file with proper headers and consistent formatting.

    Args:
        align_results: List of paths to individual alignment result files.
            All files must be in TSV format with consistent columns.
        outfile: Path where the merged results will be written.
            Directory must exist and be writable.

    Raises:
        AlignmentError: If merging fails due to file I/O or format issues.
        FileNotFoundError: If any input file doesn't exist.
    """

    # Validate inputs
    if not align_results:
        raise AlignmentError("No alignment result files provided for merging")

    # Check all input files exist
    for i, result_file in enumerate(align_results):
        if not result_file.exists():
            raise FileNotFoundError(
                f"Alignment result file {i + 1} not found: {result_file}"
            )

    logger.debug(f"Merging {len(align_results)} alignment result files")

    try:
        # Read the first file to initialize the merged dataframe
        logger.debug(f"Reading initial file: {align_results[0]}")
        merged_results = pd.read_csv(
            align_results[0],
            sep="\t",
            names=CONFIG.ALIGN_COLUMNS,
            dtype=str,  # Keep as strings to avoid type conversion issues
        )

        # Append remaining files
        for i, result_file in enumerate(align_results[1:], 1):
            logger.debug(
                "Merging file %d/%d: %s", i + 1, len(align_results), result_file
            )

            additional_results = pd.read_csv(
                result_file, sep="\t", names=CONFIG.ALIGN_COLUMNS, dtype=str
            )

            merged_results = pd.concat(
                [merged_results, additional_results], ignore_index=True, copy=False
            )

        # Ensure output directory exists
        outfile.parent.mkdir(parents=True, exist_ok=True)

        # Write merged results
        logger.debug(f"Writing merged results to: {outfile}")
        merged_results.to_csv(outfile, sep="\t", header=True, index=False)

        total_alignments = len(merged_results)
        logger.info(
            "Successfully merged %d alignment results into %s",
            total_alignments,
            outfile,
        )
        logger.debug("Merge operation completed")

    except pd.errors.EmptyDataError as e:
        raise AlignmentError(f"One or more alignment files are empty: {e}")
    except pd.errors.ParserError as e:
        raise AlignmentError(f"Error parsing alignment file format: {e}")
    except OSError as e:
        raise AlignmentError(f"File I/O error during merging: {e}")
    except Exception as e:
        raise AlignmentError(f"Unexpected error during alignment merging: {e}")


def inter_pangenome_align(
    pangenome2families_seq: Dict[str, Path],
    output: Path,
    tmpdir: Path,
    identity: float = CONFIG.DEFAULT_IDENTITY,
    coverage: float = CONFIG.DEFAULT_COVERAGE,
    cov_mode: int = CONFIG.DEFAULT_COV_MODE,
    threads: int = CONFIG.DEFAULT_THREADS,
    disable_bar: bool = False,
) -> None:
    """
    Perform inter-pangenome alignment without intra-pangenome comparisons.

    This function aligns gene families between different pangenomes while
    excluding alignments within the same pangenome. It creates MMseqs2
    databases for each pangenome and performs all pairwise comparisons.

    Args:
        pangenome2families_seq:
            Dictionary mapping pangenome names to their
            respective gene family sequence files.
        output:
            Directory where alignment results will be written.
            Will be created if it doesn't exist.
        identity:
            Minimum sequence identity threshold (0.0-1.0). Defaults to 0.8.
        coverage:
            Minimum coverage threshold (0.0-1.0). Defaults to 0.8.
        cov_mode:
            Coverage mode for MMseqs2 (0-5). Defaults to 0.
        tmpdir:
            Temporary directory for operations. If None, uses system temp.
        threads:
            Number of threads for processing. Defaults to 1.
        disable_bar:
            Whether to disable progress bars. Defaults to False.

    Raises:
        AlignmentError:
            If the alignment process fails.
        AlignmentValidationError:
            If parameters are invalid.
        FileNotFoundError:
            If sequence files don't exist.

    Notes:
        Output directory and temporary directory are supposed to be already validated.
        See the launch function to see how validation is done.
    """
    # Validate inputs
    if len(pangenome2families_seq) < 2:
        raise AlignmentValidationError(
            "At least 2 pangenomes required for inter-pangenome alignment"
        )

    logger.info(f"Processing {len(pangenome2families_seq)} pangenomes")

    pangenome2db = {}

    try:
        # Create MMseqs2 databases for each pangenome
        logger.debug("Creating MMseqs2 databases for pangenomes...")

        for name, sequences_file in pangenome2families_seq.items():
            if not sequences_file.exists():
                raise FileNotFoundError(
                    f"Sequence file not found for {name}: {sequences_file}"
                )

            logger.debug(f"Creating database for pangenome: {name}")
            try:
                # Create the database using the refactored createdb function
                pangenome2db[name] = createdb(
                    seq_files=[sequences_file], output=tmpdir, db_name=f"db_{name}"
                )
                logger.debug(f"Database created for {name}: {pangenome2db[name]}")

            except (PangenomeProcessingError, MMSeqsError) as e:
                raise AlignmentError(f"Failed to create database for {name}: {e}")

        # Perform all pairwise alignments
        logger.info("Performing pairwise alignments between pangenomes...")
        align_results = align_pangenomes(
            pangenome2db=pangenome2db,
            tmpdir=tmpdir,
            identity=identity,
            coverage=coverage,
            cov_mode=cov_mode,
            threads=threads,
            disable_bar=disable_bar,
        )

        # Merge all alignment results
        logger.debug("Merging pangenome gene families alignment results...")
        outfile = output / CONFIG.INTER_PANGENOMES_OUTPUT
        merge_aln_res(align_results, outfile)

        logger.info(
            "Inter-pangenome gene families similarities saved to: %s",
            outfile.absolute(),
        )

    except Exception as e:
        raise AlignmentError(f"Inter-pangenome alignment process failed: {e}") from e


def all_against_all_align(
    families_seq: List[Path],
    output: Path,
    identity: float = CONFIG.DEFAULT_IDENTITY,
    coverage: float = CONFIG.DEFAULT_COVERAGE,
    cov_mode: int = CONFIG.DEFAULT_COV_MODE,
    tmpdir: Optional[Path] = None,
    threads: int = CONFIG.DEFAULT_THREADS,
) -> pd.DataFrame:
    """
    Perform all-against-all alignment of gene families including intra-pangenome comparisons.

    This function combines all gene family sequences into a single database
    and performs self-alignment, capturing both inter- and intra-pangenome
    similarities.

    Args:
        families_seq:
            List of paths to gene family sequence files from all pangenomes.
        output:
            Directory where alignment results will be written.
            Will be created if it doesn't exist.
        identity:
            Minimum sequence identity threshold (0.0-1.0). Defaults to 0.8.
        coverage:
            Minimum coverage threshold (0.0-1.0). Defaults to 0.8.
        cov_mode:
            Coverage mode for MMseqs2 (0-5). Defaults to 0.
        tmpdir:
            Temporary directory for operations. If None, uses system temp.
        threads:
            Number of threads for processing. Defaults to 1.

    Returns:
        pd.DataFrame:
            DataFrame containing all alignment results with columns
            defined in CONFIG.ALIGN_COLUMNS.

    Raises:
        AlignmentError:
            If the alignment process fails.
        AlignmentValidationError:
            If parameters are invalid.
        FileNotFoundError:
            If sequence files don't exist.

    Notes:
        Output directory and temporary directory are supposed to be already validated.
        See the launch function to see how validation is done.
    """
    _validate_alignment_parameters(identity, coverage, cov_mode, threads)

    # Validate all sequence files exist
    for i, seq_file in enumerate(families_seq):
        if not seq_file.exists():
            raise FileNotFoundError(
                f"Gene family sequence file {i + 1} not found: {seq_file}"
            )

    logger.info("Starting all-against-all gene families alignment...")
    logger.info(f"Processing {len(families_seq)} sequence files")

    try:
        # Create combined MMseqs2 database
        logger.debug("Creating combined MMseqs2 database from all sequences...")
        merged_db = createdb(seq_files=families_seq, output=tmpdir)
        logger.debug(f"Combined database created: {merged_db}")

        # Create the temporary alignment database
        aln_db = tmpdir / "aln_db"
        logger.debug("Performing self-alignment on combined database...")

        # Perform self-alignment
        align_db(
            query_db=merged_db,
            target_db=merged_db,
            tmpdir=tmpdir,
            query_name="all",
            target_name="all",
            aln_db=aln_db,
            identity=identity,
            coverage=coverage,
            cov_mode=cov_mode,
            threads=threads,
        )

        # Convert alignment results to readable format
        aln_results_file = tmpdir / "all_against_all_results.tsv"

        logger.debug("Converting alignment results to TSV format...")
        write_alignment(
            query_db=merged_db,
            target_db=merged_db,
            aln_db=aln_db,
            outfile=aln_results_file,
            threads=threads,
        )

        logger.debug("Reading alignment results...")
        # Read alignment results into DataFrame
        alignment_df = pd.read_csv(
            aln_results_file,
            sep="\t",
            names=CONFIG.ALIGN_COLUMNS,
            dtype=str,  # Keep as strings initially to avoid conversion issues
        )

        # Convert numeric columns appropriately
        numeric_columns = [
            "identity",
            "qlength",
            "tlength",
            "alnlength",
            "e_value",
            "bits",
        ]
        for col in numeric_columns:
            if col in alignment_df.columns:
                alignment_df[col] = pd.to_numeric(alignment_df[col], errors="coerce")

        # Save results to the output file
        outfile = output / CONFIG.ALL_AGAINST_ALL_OUTPUT
        logger.debug(f"Saving results to: {outfile}")
        alignment_df.to_csv(outfile, sep="\t", header=True, index=False)

        total_alignments = len(alignment_df)
        logger.info(
            "All-against-all alignment completed with %d results saved to: %s",
            total_alignments,
            outfile.absolute(),
        )

        return alignment_df

    except (PangenomeProcessingError, MMSeqsError) as e:
        raise AlignmentError(f"Database creation failed: {e}") from e
    except pd.errors.EmptyDataError as e:
        raise AlignmentError("Alignment produced no results") from e
    except Exception as e:
        raise AlignmentError(f"All-against-all alignment process failed: {e}") from e


def launch_pangenomes_alignment(
    pangenomes: Pangenomes,
    output: Path,
    tmpdir: Path,
    inter_pangenomes: bool = False,
    all_against_all: bool = False,
    identity: float = CONFIG.DEFAULT_IDENTITY,
    coverage: float = CONFIG.DEFAULT_COVERAGE,
    cov_mode: int = CONFIG.DEFAULT_COV_MODE,
    threads: int = CONFIG.DEFAULT_THREADS,
    lock: Lock = None,
    force: bool = False,
    disable_bar: bool = False,
) -> None:
    """
    Launches the alignment of pangenome families based on the specified mode.

    Args:
        pangenomes (Pangenomes):
            A collection of pangenomes to be aligned.
        output (Path):
            Path to the directory where the alignment results will be stored.
        tmpdir (Path):
            Path to the temporary directory for intermediate files during
            the alignment process.
        inter_pangenomes (bool):
            If True, performs an inter-pangenome alignment mode.
        all_against_all (bool):
            If True, performs an all-against-all alignment mode.
        identity (float):
            Minimum sequence identity threshold for the alignment.
        coverage (float):
            Minimum sequence coverage threshold for the alignment.
        cov_mode (int):
            Coverage mode to dictate how the coverage threshold is applied.
        threads (int):
            Number of threads to be used for parallel processing.
        lock (Lock):
            A multiprocessing lock to synchronize access to certain operations.
        force (bool):
            If True, allows overwriting or recreating the output directory.
        disable_bar (bool):
            If True, disables progress bars in the alignment process.

    Raises:
        AlignmentValidationError:
            If none of the alignment modes (-inter_pangenomes
            or -all_against_all) are specified.
    """
    # Write pangenome families sequences
    logger.info("Writing pangenome families sequences...")
    pangenome2families_seq = write_pangenomes_families_sequences(
        pangenomes=pangenomes,
        output=tmpdir,
        threads=threads,
        lock=init_lock(lock),
        disable_bar=disable_bar,
    )

    logger.info(
        f"Successfully wrote sequences for {len(pangenome2families_seq)} pangenomes"
    )

    # Create output directory
    logger.debug(f"Creating output directory: {output}")
    mkdir(output, force=force)

    # Perform alignment based on the selected mode
    if inter_pangenomes:
        logger.info("Starting inter-pangenome alignment mode...")
        inter_pangenome_align(
            pangenome2families_seq=pangenome2families_seq,
            output=output,
            identity=identity,
            coverage=coverage,
            cov_mode=cov_mode,
            tmpdir=tmpdir,
            threads=threads,
            disable_bar=disable_bar,
        )

    elif all_against_all:
        logger.info("Starting all-against-all alignment mode...")
        all_against_all_align(
            families_seq=list(pangenome2families_seq.values()),
            output=output,
            identity=identity,
            coverage=coverage,
            cov_mode=cov_mode,
            tmpdir=tmpdir,
            threads=threads,
        )
    else:
        raise AlignmentValidationError(
            "You must choose between --inter_pangenomes or --all_against_all alignment modes"
        )


def launch(args: argparse.Namespace) -> None:
    """
    Main entry point for alignment operations.

    This function orchestrates the complete alignment workflow, including
    parameter validation, pangenome loading, sequence writing, and alignment
    execution based on the specified mode (inter-pangenome or all-against-all).

    Args:
        args: Parsed command line arguments containing all configuration
            parameters for the alignment operation.

    Raises:
        AlignmentError: If any step of the alignment process fails.
        AlignmentValidationError: If parameter validation fails.
        NotADirectoryError: If specified directories are invalid.
    """
    # Validate command line parameters
    logger.info("Validating alignment parameters...")
    check_align_parameters(args)

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
        check_function=check_pangenome_align,
        max_workers=args.threads,
        lock=lock,
        disable_bar=args.disable_prog_bar,
    )

    logger.info(f"Successfully loaded {len(pangenomes)} pangenomes")

    # Create a temporary directory
    tmpdir = Path(tempfile.mkdtemp(dir=args.tmpdir))
    logger.debug(f"Created temporary directory: {tmpdir}")

    try:
        launch_pangenomes_alignment(
            pangenomes=pangenomes,
            output=args.output,
            tmpdir=tmpdir,
            inter_pangenomes=args.inter_pangenomes,
            all_against_all=args.all_against_all,
            identity=args.align_identity,
            coverage=args.align_coverage,
            cov_mode=args.align_cov_mode,
            threads=args.threads,
            lock=lock,
            force=args.force,
            disable_bar=args.disable_prog_bar,
        )
        logger.info("Alignment process completed successfully")
    except Exception as e:
        raise AlignmentError(f"Alignment process failed: {e}") from e
    finally:
        # Clean up the temporary directory unless requested to keep it
        if not args.keep_tmp:
            logger.debug(f"Cleaning up temporary directory: {tmpdir}")
            rmtree(tmpdir, ignore_errors=True)
        else:
            logger.info(f"Temporary files kept in: {tmpdir}")


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Create argument subparser for alignment command.

    This function sets up the command-line interface for the alignment
    functionality within the PANORAMA tool suite.

    Args:
        sub_parser: The subparser object from argparse to add alignment command to.

    Returns:
        argparse.ArgumentParser: Configured parser for alignment command.
    """
    parser = sub_parser.add_parser(
        "align",
        help="Align gene families between pangenomes using MMseqs2",
        description="Perform sequence alignment between pangenome gene families "
        "using MMseqs2 with support for both inter-pangenome and "
        "all-against-all alignment modes.",
    )
    parser_align(parser)
    return parser


def parser_mmseqs2_align(parser: argparse.ArgumentParser) -> argparse._ArgumentGroup:
    """
    Add MMseqs2-specific arguments to the parser.

    Args:
        parser: ArgumentParser to add MMseqs2 arguments to.

    Returns:
        argparse._ArgumentGroup: The argument group containing MMseqs2 options.
    """
    mmseqs_group = parser.add_argument_group(
        title="MMseqs2 alignment parameters",
        description="Configure MMseqs2 alignment behavior. "
        "See MMseqs2 documentation for detailed parameter descriptions.",
    )

    mmseqs_group.add_argument(
        "--align_identity",
        required=False,
        type=float,
        default=CONFIG.DEFAULT_IDENTITY,
        metavar="FLOAT",
        help=f"Minimum identity percentage threshold (0.0-1.0). Default: {CONFIG.DEFAULT_IDENTITY}",
    )

    mmseqs_group.add_argument(
        "--align_coverage",
        required=False,
        type=float,
        default=CONFIG.DEFAULT_COVERAGE,
        metavar="FLOAT",
        help=f"Minimum coverage percentage threshold (0.0-1.0). Default: {CONFIG.DEFAULT_COVERAGE}",
    )

    mmseqs_group.add_argument(
        "--align_cov_mode",
        required=False,
        type=int,
        default=CONFIG.DEFAULT_COV_MODE,
        choices=[0, 1, 2, 3, 4, 5],
        metavar="INT",
        help=f"Coverage mode: 0=query, 1=target, 2=shorter seq, 3=longer seq, "
        f"4=query and target, 5=shorter and longer seq. Default: {CONFIG.DEFAULT_COV_MODE}",
    )

    return mmseqs_group


def parser_align(parser: argparse.ArgumentParser) -> None:
    """
    Configure the argument parser for the alignment command.

    This function adds all necessary command-line arguments for the alignment
    functionality, including required arguments, alignment modes, and optional
    parameters.

    Args:
        parser: ArgumentParser to configure with alignment arguments.
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
        help="Output directory where alignment results will be written",
    )

    # Mutually exclusive alignment mode selection
    exclusive_group = required_group.add_mutually_exclusive_group(required=True)

    exclusive_group.add_argument(
        "--inter_pangenomes",
        action="store_true",
        help="Perform inter-pangenome alignment only (exclude intra-pangenome comparisons). "
        "Cannot be used with --all_against_all",
    )

    exclusive_group.add_argument(
        "--all_against_all",
        action="store_true",
        help="Perform all-against-all alignment including intra-pangenome comparisons. "
        "Cannot be used with --inter_pangenomes",
    )

    # Add MMseqs2 specific arguments
    parser_mmseqs2_align(parser)

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
