#!/usr/bin/env python3
# coding:utf-8

"""
Module for creating MMseqs2 databases and processing pangenome families sequences.

This module provides functionality to create MMseqs2 sequence databases and write
pangenome families sequences using multithreading for improved performance.
"""

# default libraries
from __future__ import annotations
import logging
import tempfile
import subprocess

from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Lock
from dataclasses import dataclass
from contextlib import contextmanager

# installed libraries
from tqdm import tqdm
from ppanggolin.formats.writeSequences import write_fasta_prot_fam_from_pangenome_file

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.utils import init_lock, mkdir


logger = logging.getLogger("PANORAMA")


# Configuration constants
@dataclass(frozen=True)
class MMSeqsConfig:
    """Configuration constants for MMseqs2 operations."""

    DEFAULT_DB_TYPE: int = 0
    PROTEIN_FAMILIES_FILENAME: str = "all_protein_families.faa.gz"
    FAMILY_FILTER_ALL: str = "all"


class PangenomeProcessingError(Exception):
    """Custom exception for pangenome processing errors."""

    pass


class MMSeqsError(Exception):
    """Custom exception for MMseqs2 related errors."""

    pass


def _validate_seq_files(seq_files: List[Path]) -> None:
    """
    Validate sequence files exist and are readable.

    Args:
        seq_files: List of sequence file paths to validate.

    Raises:
        PangenomeProcessingError: If any file is invalid or doesn't exist.
    """
    if not seq_files:
        raise PangenomeProcessingError("No sequence files provided")

    for seq_file in seq_files:
        if not isinstance(seq_file, Path):
            raise PangenomeProcessingError(f"Invalid file type: {type(seq_file)}")
        if not seq_file.exists():
            raise PangenomeProcessingError(f"Sequence file does not exist: {seq_file}")
        if not seq_file.is_file():
            raise PangenomeProcessingError(f"Path is not a file: {seq_file}")


def createdb(
    seq_files: List[Path],
    output: Path,
    db_name: str = "mmseqs_db",
    db_type: int = MMSeqsConfig.DEFAULT_DB_TYPE,
) -> Path:
    """
    Create a MMseqs2 sequence database from the given FASTA files.

    This function creates an MMseqs2 database by combining multiple sequence files
    into a single database that can be used for sequence searches and clustering.

    Args:
        seq_files: List of FASTA file paths to include in the database.
            All files must exist and be readable.
        output: Temporary directory where the database will be created.
            Directory must exist and be writable.
        db_type: Type of MMseqs2 database to create. Defaults to 0 (protein).
            Valid values: 0 (protein), 1 (nucleotide), 2 (HMM profile).

    Returns:
        Path: Path to the created MMseqs2 database file.

    Raises:
        PangenomeProcessingError: If input validation fails or file operations fail.
        MMSeqsError: If MMseqs2 command execution fails.
    """
    # Input validation
    _validate_seq_files(seq_files)
    mkdir(output, force=True)

    if not isinstance(db_type, int) or db_type < 0:
        raise PangenomeProcessingError(
            f"Invalid db_type: {db_type}. Must be non-negative integer."
        )

    logger.debug(
        f"Creating MMseqs2 database from {len(seq_files)} sequence files",
    )

    # Create a database using a context manager for proper resource management
    db_path = Path(output / f"{db_name}")

    try:
        # Prepare MMseqs2 command
        absolute_paths = [str(seq_file.absolute()) for seq_file in seq_files]
        cmd = [
            "mmseqs",
            "createdb",
            *absolute_paths,
            db_path.absolute().as_posix(),
            "--dbtype",
            str(db_type),
        ]

        logger.debug(f"Executing MMseqs2 command: {' '.join(cmd)}")

        # Execute MMseqs2 command with proper error handling
        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=False,
            text=True,
        )

        if result.returncode != 0:
            error_msg = f"MMseqs2 createdb failed with return code {result.returncode}"
            if result.stderr:
                error_msg += f": {result.stderr.strip()}"
            raise MMSeqsError(error_msg)

        logger.debug(f"MMseqs2 database successfully created at: {db_path}")
        return db_path

    except subprocess.SubprocessError as e:
        raise MMSeqsError(f"Failed to execute MMseqs2 command: {e}")
    except Exception as e:
        raise PangenomeProcessingError(
            f"Unexpected error during database creation: {e}"
        )


def write_protein_families_sequences(
    pangenome: Pangenome, output_dir: Path
) -> Tuple[str, Path]:
    """
    Write protein families sequences for a single pangenome.

    This is a wrapper function designed for use in multithreading contexts.

    Args:
        pangenome: Pangenome object containing gene families data.
        output_dir: Directory where sequences will be written.

    Returns:
        Tuple[str, Path]: Pangenome name and path to the generated sequences file.

    Raises:
        PangenomeProcessingError: If sequence writing fails.

    Note:
        This function prefixes sequence names with pangenome name to ensure
        uniqueness across different pangenomes when combining sequences.
    """
    try:
        # Configure output parameters
        write_params = {
            "output": output_dir,
            "family_filter": MMSeqsConfig.FAMILY_FILTER_ALL,
            "compress": True,
            "disable_bar": True,
        }

        # Write sequences using the ppanggolin function
        write_fasta_prot_fam_from_pangenome_file(pangenome.file, **write_params)

        # Return pangenome name and path to the generated file
        sequences_path = output_dir / MMSeqsConfig.PROTEIN_FAMILIES_FILENAME

        if not sequences_path.exists():
            raise PangenomeProcessingError(
                f"Expected output file not created: {sequences_path}"
            )

        return pangenome.name, sequences_path

    except Exception as e:
        raise PangenomeProcessingError(
            f"Failed to write sequences for pangenome {pangenome.name}: {e}"
        )


def write_pangenomes_families_sequences(
    pangenomes: Pangenomes,
    output: Path,
    threads: int = 1,
    lock: Optional[Lock] = None,
    disable_bar: bool = False,
) -> Dict[str, Path]:
    """
    Write protein families sequences for multiple pangenomes using multithreading.

    This function processes multiple pangenomes concurrently to extract and write
    protein families sequences. Each pangenome's sequences are written to a separate
    subdirectory within the temporary directory.

    Args:
        pangenomes: Pangenomes object containing multiple pangenome instances.
            Must contain at least one pangenome.
        output: Output directory where sequence files will be written.
            Directory will be created if it doesn't exist.
        threads: Number of worker threads to use for parallel processing.
            Defaults to 1. Must be a positive integer.
        lock: Optional multiprocessing Lock object for thread-safe operations.
            If None, a new lock will be initialized for thread safety.
        disable_bar: Whether to disable the progress bar display.
            Defaults to False (progress bar will be shown).

    Returns:
        Dict[str, Path]: Dictionary mapping pangenome names to their corresponding
            protein families sequence file paths.

    Raises:
        PangenomeProcessingError: If input validation fails or processing errors occur.

    Note:
        - Each pangenome's sequences are written to tmpdir/{pangenome_name}/
        - Sequence names are prefixed with pangenome name to avoid duplicates
        - Files are compressed (.gz) to save space
    """
    # Input validation
    mkdir(output, force=True, erase=True)

    logger.info(
        f"Writing families sequences for {len(pangenomes)} pangenomes using {threads} threads..."
    )

    # Initialize results dictionary
    pangenomes2families_sequences: Dict[str, Path] = {}

    try:
        # Use ThreadPoolExecutor with proper initialization
        with ThreadPoolExecutor(
            max_workers=threads, initializer=init_lock, initargs=(lock,)
        ) as executor:

            # Submit all tasks
            future_to_pangenome = {}
            for pangenome in pangenomes:
                # Create an output directory for this pangenome
                output_dir = mkdir(output / f"{pangenome.name}")

                # Submit task to executor
                future = executor.submit(
                    write_protein_families_sequences, pangenome, output_dir
                )
                future_to_pangenome[future] = pangenome.name

            # Process completed tasks with progress tracking
            with tqdm(
                total=len(pangenomes),
                unit="pangenome",
                disable=disable_bar,
                desc="Processing pangenomes",
            ) as progress_bar:

                # Use as_completed for better error handling and progress tracking
                for future in as_completed(future_to_pangenome):
                    pangenome_name = future_to_pangenome[future]

                    try:
                        # Get result from a completed task
                        name, sequences_path = future.result()
                        pangenomes2families_sequences[name] = sequences_path

                        logger.debug(f"Completed processing pangenome: {name}")
                    except Exception as e:
                        raise PangenomeProcessingError(
                            f"Error processing pangenome {pangenome_name}: {e}"
                        )
                    finally:
                        # Update the progress bar regardless of success/failure
                        progress_bar.update(1)
        logger.info(
            f"Successfully processed {len(pangenomes2families_sequences)} pangenomes"
        )
        return pangenomes2families_sequences

    except Exception as e:
        raise PangenomeProcessingError(
            f"Failed to write pangenomes families sequences: {e}"
        )
