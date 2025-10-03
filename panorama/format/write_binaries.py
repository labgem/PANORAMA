#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to write, update and erase pangenome system data from HDF5 files.

The module extends the base ppanggolin functionality to handle system-specific data structures
including canonical systems, system units and their metadata relationships.
"""

# default libraries
import logging
from dataclasses import dataclass
from typing import Dict, Set, Tuple, Union

# Third-party imports
import tables

# PPanGGOLiN format imports
from ppanggolin.formats.writeBinaries import erase_pangenome as super_erase_pangenome
from ppanggolin.formats.writeBinaries import write_pangenome as super_write_pangenome
from ppanggolin.formats.writeBinaries import write_status as super_write_status
from tqdm import tqdm

# local libraries
from panorama.pangenomes import Pangenome

# Module-level constants
DEFAULT_STRING_SIZE = 1
DEFAULT_INT_SIZE = tables.Int64Col()

# Create module logger
logger = logging.getLogger("PANORAMA")


@dataclass
class SystemTableSizes:
    """Container for system table size information."""

    system_sizes: Tuple[int, int, int]  # (max_id_len, max_name_len, expected_rows)
    unit_sizes: Tuple[int, int, int, int]  # (max_name_len, max_gf_name_len, max_metadata_source, expected_rows)
    canonical_sizes: Tuple[int, int, int]  # (max_id_len, max_name_len, expected_rows)
    canonical_unit_sizes: Tuple[
        int, int, int, int
    ]  # (max_name_len, max_gf_name_len, max_metadata_source, expected_rows)
    system_to_canonical_rows: int


def create_system_unit_description(
    max_name_len: int = DEFAULT_STRING_SIZE,
    max_gf_name_len: int = DEFAULT_STRING_SIZE,
    max_metadata_source: int = DEFAULT_STRING_SIZE,
) -> Dict[str, Union[tables.StringCol, tables.Int64Col]]:
    """
    Create the dictionary that describes the HDF5 table structure for system units.

    Args:
        max_name_len (int): Maximum size of system unit name. Defaults to 1.
        max_gf_name_len (int): Maximum size of gene family name. Defaults to 1.
        max_metadata_source (int): Maximum size of annotation source name. Defaults to 1.

    Returns:
        Dict[str, Union[tables.StringCol, tables.Int64Col]]: HDF5 table description dictionary
            containing column definitions for the system unit table.
    """
    return {
        "ID": DEFAULT_INT_SIZE,
        "name": tables.StringCol(itemsize=max_name_len),
        "geneFam": tables.StringCol(itemsize=max_gf_name_len),
        "metadata_id": DEFAULT_INT_SIZE,
        "metadata_source": tables.StringCol(itemsize=max_metadata_source),
    }


def create_system_description(
    max_id_len: int = DEFAULT_STRING_SIZE,
    max_name_len: int = DEFAULT_STRING_SIZE,
    include_canonical_count: bool = False,
) -> Dict[str, Union[tables.StringCol, tables.Int64Col]]:
    """
    Create the dictionary that describes the HDF5 table structure for systems.

    Args:
        max_id_len (int): Maximum size of system identifier. Defaults to 1.
        max_name_len (int): Maximum size of system name. Defaults to 1.
        include_canonical_count (bool): Whether to include the canonical count column. Defaults to False.

    Returns:
        Dict[str, Union[tables.StringCol, tables.Int64Col]]: HDF5 table description dictionary
            containing column definitions for the system table.
    """
    description = {
        "ID": tables.StringCol(itemsize=max_id_len),
        "name": tables.StringCol(itemsize=max_name_len),
        "unit": DEFAULT_INT_SIZE,
    }

    # Add the canonical count column if requested
    if include_canonical_count:
        description["canonical"] = DEFAULT_INT_SIZE

    return description


def _calculate_max_lengths_for_system(
    system, current_max_sys: Tuple[int, int], current_max_unit: Tuple[int, int, int]
) -> Tuple[Tuple[int, int], Tuple[int, int, int]]:
    """
    Calculate maximum string lengths for a single system and its units.

    Args:
        system: System object to analyze
        current_max_sys (Tuple[int, int]): Current maximum (id_len, name_len) for systems
        current_max_unit (Tuple[int, int, int]): Current maximum (name_len, gf_name_len, annot_source_len) for units

    Returns:
        Tuple[Tuple[int, int], Tuple[int, int, int]]: Updated maximum lengths for systems and units
    """
    max_id_len, max_sys_name_len = current_max_sys
    max_unit_name_len, max_gf_name_len, max_annot_source_len = current_max_unit

    # Update system-level maximums
    max_id_len = max(max_id_len, len(system.ID))
    max_sys_name_len = max(max_sys_name_len, len(system.name))

    # Update annotation source maximums
    for annot_source in system.annotation_sources():
        max_annot_source_len = max(max_annot_source_len, len(annot_source))

    # Update unit-level maximums
    for unit in system.units:
        max_unit_name_len = max(max_unit_name_len, len(unit.name))
        for gene_family in unit.families:
            max_gf_name_len = max(max_gf_name_len, len(gene_family.name))

    return (max_id_len, max_sys_name_len), (
        max_unit_name_len,
        max_gf_name_len,
        max_annot_source_len,
    )


def calculate_system_table_sizes(pangenome: Pangenome, source: str) -> SystemTableSizes:
    """
    Calculate maximum sizes and expected row counts for all system-related HDF5 tables.

    This function analyzes all systems from a given source to determine the optimal
    table sizes for efficient HDF5 storage.

    Args:
        pangenome (Pangenome): Pangenome object containing the systems data
        source (str): Name of the system source to analyze

    Returns:
        SystemTableSizes: Dataclass containing all calculated size information for
            system tables, unit tables, canonical tables and cross-reference tables.
    """
    # Initialize tracking variables
    max_len_sys = (DEFAULT_STRING_SIZE, DEFAULT_STRING_SIZE)
    max_len_unit = (DEFAULT_STRING_SIZE, DEFAULT_STRING_SIZE, DEFAULT_STRING_SIZE)
    max_len_canonical = (DEFAULT_STRING_SIZE, DEFAULT_STRING_SIZE)
    max_len_canon_unit = (DEFAULT_STRING_SIZE, DEFAULT_STRING_SIZE, DEFAULT_STRING_SIZE)

    # Initialize row counters
    exp_rows_sys = exp_rows_units = exp_rows_can = exp_rows_can_units = exp_rows_sys2can = 0

    # Analyze each system to calculate maximums and counts
    for system in pangenome.get_system_by_source(source):
        # Update maximum lengths for regular systems
        max_len_sys, max_len_unit = _calculate_max_lengths_for_system(system, max_len_sys, max_len_unit)

        # Count rows for regular systems
        exp_rows_sys += len(system)
        exp_rows_units += system.number_of_families

        # Process canonical systems
        for canonical in system.canonical:
            max_len_canonical, max_len_canon_unit = _calculate_max_lengths_for_system(
                canonical, max_len_canonical, max_len_canon_unit
            )
            exp_rows_can += len(canonical)
            exp_rows_can_units += canonical.number_of_families
            exp_rows_sys2can += 1  # One row per system-canonical relationship

    return SystemTableSizes(
        system_sizes=(*max_len_sys, exp_rows_sys),
        unit_sizes=(*max_len_unit, exp_rows_units),
        canonical_sizes=(*max_len_canonical, exp_rows_can),
        canonical_unit_sizes=(*max_len_canon_unit, exp_rows_can_units),
        system_to_canonical_rows=exp_rows_sys2can,
    )


def write_system_data_to_tables(
    system,
    system_row: tables.Table.row,
    unit_row: tables.Table.row,
    is_canonical: bool = False,
) -> None:
    """
    Write a single system's data to the appropriate HDF5 table rows.

    Args:
        system: System object containing the data to write
        system_row (tables.Table.row): HDF5 table row object for system data
        unit_row (tables.Table.row): HDF5 table row object for unit data
        is_canonical (bool): Whether this system is a canonical version. Defaults to False.

    Returns:
        None

    Note:
        This function modifies the table rows in-place and calls append() to commit the data.
    """
    # Write each unit of the system
    for unit in system.units:
        # Populate system-level data
        system_row["ID"] = system.ID
        system_row["name"] = system.name
        system_row["unit"] = unit.ID

        # Add canonical count only for non-canonical systems
        if not is_canonical:
            system_row["canonical"] = len(system.canonical)

        system_row.append()

        # Write unit-level data for each gene family in the system
        for gene_family in system.families:
            unit_row["ID"] = unit.ID
            unit_row["name"] = unit.name
            unit_row["geneFam"] = gene_family.name

            # Extract metadata information
            metadata_source, metadata_id = unit.get_metainfo(gene_family)
            unit_row["metadata_source"] = metadata_source
            unit_row["metadata_id"] = metadata_id

            unit_row.append()


def create_system_tables(
    h5f: tables.File, source_group: tables.Group, table_sizes: SystemTableSizes
) -> Tuple[tables.Table, ...]:
    """
    Create all HDF5 tables needed for storing system data.

    Args:
        h5f (tables.File): HDF5 file handles
        source_group (tables.Group): HDF5 group for this source
        table_sizes (SystemTableSizes): Size information for table creation

    Returns:
        Tuple[tables.Table, ...]: Tuple of created tables (system, unit, canonical, canonical_unit, sys2canonical)
    """
    # Create the main system table
    system_table = h5f.create_table(
        source_group,
        "systems",
        description=create_system_description(*table_sizes.system_sizes[:-1], include_canonical_count=True),
        expectedrows=table_sizes.system_sizes[-1],
    )

    # Create the unit table
    unit_table = h5f.create_table(
        source_group,
        "units",
        description=create_system_unit_description(*table_sizes.unit_sizes[:-1]),
        expectedrows=table_sizes.unit_sizes[-1],
    )

    # Create the canonical system table
    canonical_table = h5f.create_table(
        source_group,
        "canonic",
        description=create_system_description(*table_sizes.canonical_sizes[:-1]),
        expectedrows=table_sizes.canonical_sizes[-1],
    )

    # Create the canonical unit table
    canonical_unit_table = h5f.create_table(
        source_group,
        "canonic_units",
        description=create_system_unit_description(*table_sizes.canonical_unit_sizes[:-1]),
        expectedrows=table_sizes.canonical_unit_sizes[-1],
    )

    # Create the system-to-canonical mapping table
    sys2canonical_table = h5f.create_table(
        source_group,
        "system_to_canonical",
        description={
            "system": tables.StringCol(itemsize=table_sizes.system_sizes[0]),
            "canonic": tables.StringCol(itemsize=table_sizes.canonical_sizes[0]),
        },
        expectedrows=table_sizes.system_to_canonical_rows,
    )

    return (
        system_table,
        unit_table,
        canonical_table,
        canonical_unit_table,
        sys2canonical_table,
    )


def write_systems_to_hdf5(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False) -> None:
    """
    Write all systems from a specific source to HDF5 file.

    This function creates the necessary HDF5 table structure and writes all system data,
    including regular systems, canonical systems and their relationships.

    Args:
        pangenome (Pangenome): Pangenome object containing systems to write
        h5f (tables.File): Open HDF5 file handle for writing
        source (str): Source identifier for the systems being written
        disable_bar (bool): Whether to disable the progress bar. Defaults to False.

    Returns:
        None

    Raises:
        AssertionError: If systems data is not properly formatted or accessible
    """
    # Ensure systems group exists
    if "/systems" not in h5f:
        systems_group = h5f.create_group("/", "systems", "Detected systems")
    else:
        systems_group = h5f.root.systems

    # Calculate optimal table sizes
    logger.debug(f"Calculating table sizes for systems from source: {source}")
    table_sizes = calculate_system_table_sizes(pangenome, source)

    # Create the source-specific group
    source_group = h5f.create_group(systems_group, source, f"Detected systems from source: {source}")

    # Store metadata sources information
    source_group._v_attrs.metadata_sources = pangenome.systems_sources_to_metadata_source()[source]

    # Create all necessary tables
    tables_tuple = create_system_tables(h5f, source_group, table_sizes)
    (
        system_table,
        unit_table,
        canonical_table,
        canonical_unit_table,
        sys2canonical_table,
    ) = tables_tuple

    # Get table row objects for efficient writing
    system_row = system_table.row
    unit_row = unit_table.row
    canonical_row = canonical_table.row
    canonical_unit_row = canonical_unit_table.row
    sys2canonical_row = sys2canonical_table.row

    # Track canonical systems to avoid duplicates
    canonical_systems_seen: Set[str] = set()

    # Write all systems with progress tracking
    total_systems = pangenome.number_of_systems(source=source, with_canonical=False)
    with tqdm(
        total=total_systems,
        unit="system",
        disable=disable_bar,
        desc=f"Writing systems from {source}",
    ) as progress:
        for system in pangenome.systems:
            # Write regular system data
            write_system_data_to_tables(system, system_row, unit_row, is_canonical=False)

            # Process canonical systems
            for canonical in system.canonical:
                # Write canonical system data only once
                if canonical.ID not in canonical_systems_seen:
                    canonical_systems_seen.add(canonical.ID)
                    write_system_data_to_tables(canonical, canonical_row, canonical_unit_row, is_canonical=True)

                # Write system-to-canonical relationship
                sys2canonical_row["system"] = system.ID
                sys2canonical_row["canonic"] = canonical.ID
                sys2canonical_row.append()

            progress.update()

    # Flush all tables to ensure data is written
    for table in [
        system_table,
        unit_table,
        canonical_table,
        canonical_unit_table,
        sys2canonical_table,
    ]:
        table.flush()

    logger.debug(f"Successfully wrote {total_systems} systems with {len(canonical_systems_seen)} canonical systems")


def write_pangenome_status(pangenome: Pangenome, h5f: tables.File) -> None:
    """
    Write pangenome processing status to the HDF5 file.

    This function extends the base status writing functionality to include
    system-specific status information.

    Args:
        pangenome (Pangenome): Pangenome object with current status
        h5f (tables.File): Open HDF5 file handle

    Returns:
        None
    """
    # Call parent status writing function
    super_write_status(pangenome, h5f)

    status_group = h5f.root.status

    # Handle systems status
    if "systems" in pangenome.status:
        # Initialize systems_sources attribute if not present
        if not hasattr(status_group._v_attrs, "systems_sources"):
            status_group._v_attrs.systems_sources = set()

        # Update status based on the current state
        if pangenome.status["systems"] in ["Computed", "Loaded"]:
            status_group._v_attrs.systems = True
            status_group._v_attrs.systems_sources |= pangenome.systems_sources
        else:
            status_group._v_attrs.systems = False
    else:
        # No systems data present
        status_group._v_attrs.systems = False
        status_group._v_attrs.systems_sources = set()


def erase_pangenome(
    pangenome: Pangenome,
    graph: bool = False,
    gene_families: bool = False,
    partition: bool = False,
    rgp: bool = False,
    spots: bool = False,
    modules: bool = False,
    metadata: bool = False,
    systems: bool = False,
    source: str = None,
) -> None:
    """
    Erase specific data tables from a pangenome HDF5 file.

    This function provides selective deletion of pangenome data components,
    extending the base functionality to handle system-specific data.

    Args:
        pangenome (Pangenome): Pangenome object to modify
        graph (bool): Remove graph information. Defaults to False.
        gene_families (bool): Remove gene families information. Defaults to False.
        partition (bool): Remove partition information. Defaults to False.
        rgp (bool): Remove RGP information. Defaults to False.
        spots (bool): Remove spots information. Defaults to False.
        modules (bool): Remove modules' information. Defaults to False.
        metadata (bool): Remove metadata. Defaults to False.
        systems (bool): Remove systems data. Defaults to False.
        source (str): Specific source to remove (required for systems/metadata). Defaults to None.

    Returns:
        None

    Raises:
        AssertionError: If the source is None when systems=True
        FileNotFoundError: If the pangenome file doesn't exist
    """
    # Open the HDF5 file for modification
    with tables.open_file(pangenome.file, "a") as h5f:
        status_group = h5f.root.status

        # Call parent erase function for standard components
        super_erase_pangenome(
            pangenome,
            graph,
            gene_families,
            partition,
            rgp,
            spots,
            modules,
            metadata=metadata,
            metatype="families" if metadata else None,
            source=source,
        )

        # Handle systems-specific erasure
        if "/systems" in h5f and systems:
            if source is None:
                raise ValueError("Source must be specified when erasing systems data")

            systems_group = h5f.root.systems

            # Remove the specific source if it exists
            if source in systems_group:
                logger.info(f"Erasing systems data from source: {source}")
                h5f.remove_node("/systems", source, recursive=True)

                # Update status tracking
                status_group._v_attrs.systems_sources.discard(source)
                pangenome.status["systems_sources"].discard(source)

                # Remove the entire systems group if no sources remain
                if len(status_group._v_attrs.systems_sources) == 0:
                    h5f.remove_node("/", "systems")
                    status_group._v_attrs.systems = False
                    pangenome.status["systems"] = "No"
            else:
                logger.warning(f"Systems source '{source}' not found in HDF5 file")


def write_pangenome(
    pangenome: Pangenome,
    file_path: str,
    source: str = None,
    force: bool = False,
    disable_bar: bool = False,
) -> None:
    """
    Write or update a complete pangenome to an HDF5 file.

    This function handles the complete workflow of writing pangenome data,
    including both standard components and system-specific data.

    Args:
        pangenome (Pangenome): Pangenome object containing all data to write
        file_path (str): Path to the HDF5 file for the output
        source (str): Source identifier for systems or metadata. Defaults to None.
        force (bool): Whether to overwrite existing files. Defaults to False.
        disable_bar (bool): Whether to disable progress bars. Defaults to False.

    Returns:
        None

    Raises:
        AssertionError: If the source is None when systems need to be written
        IOError: If the file cannot be created or accessed
    """
    # Write standard pangenome components first
    super_write_pangenome(pangenome, file_path, force, disable_bar)

    # Handle systems data if present and computed
    with tables.open_file(file_path, "a") as h5f:
        if "systems" in pangenome.status and pangenome.status["systems"] == "Computed":
            if source is None:
                raise ValueError("Source must be specified when writing systems data")

            logger.info(f"Writing systems data from source: {source}")
            write_systems_to_hdf5(pangenome=pangenome, h5f=h5f, source=source, disable_bar=disable_bar)

            # Update status to reflect successful writing
            pangenome.status["systems"] = "Loaded"

        # Write updated status information
        write_pangenome_status(pangenome, h5f)

    logger.info(f"Successfully wrote pangenome to: {file_path}")
