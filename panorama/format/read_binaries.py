#!/usr/bin/env python3
# coding:utf-8

"""
Pangenome Data Reader Module

This module provides comprehensive functions to read and load pangenome data from HDF5 files.
It supports parallel loading of multiple pangenomes with various data components, including
annotations, gene families, systems, spots, modules, and metadata.

The module is designed to work with the PPanGGOLiN pangenome analysis toolkit and provides
robust error handling and progress tracking for large-scale pangenome analyses.
"""

# Standard library imports
import logging
import time
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set

# Third-party imports
import tables
from tqdm import tqdm

# PPanGGOLiN format imports
from ppanggolin.formats import (
    get_need_info,
    read_annotation,
    read_chunks,
    read_gene_sequences,
    read_graph,
    read_metadata,
    read_rgp,
)
from ppanggolin.formats import get_status as super_get_status
from ppanggolin.geneFamily import Gene

# Local imports
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.region import Module, Spot
from panorama.systems.models import Models
from panorama.systems.system import System, SystemUnit
from panorama.utils import check_tsv_sanity, init_lock

# Module-level logger
logger = logging.getLogger("PANORAMA")


def get_status(pangenome: Pangenome, pangenome_file: Path) -> None:
    """
    Check which elements are present in the HDF5 file and update pangenome status.

    This function extends the base status checking functionality to include
    systems-specific status information.

    Args:
        pangenome (Pangenome): The pangenome object to update with status information.
        pangenome_file (Path): Path to the pangenome HDF5 file to examine.

    Note:
        This function modifies the pangenome object in-place by updating its status dictionary.
    """
    # Get standard status information from the parent function
    super_get_status(pangenome, pangenome_file)

    # Open the HDF5 file in read-only mode for status checking
    with tables.open_file(pangenome_file.absolute().as_posix(), "r") as h5f:
        status_group = h5f.root.status

        # Check for systems-specific status attributes
        if hasattr(status_group._v_attrs, "systems") and status_group._v_attrs.systems:
            pangenome.status["systems"] = "inFile"

            # Load systems sources if available, otherwise initialize an empty set
            if hasattr(status_group._v_attrs, 'systems_sources'):
                pangenome.status["systems_sources"] = status_group._v_attrs.systems_sources
            else:
                pangenome.status["systems_sources"] = set()


def read_systems_by_source(
        pangenome: Pangenome,
        source_group: tables.Group,
        models: Models,
        read_canonical: bool = True,
        disable_bar: bool = False
) -> None:
    """
    Read systems from a specific source and integrate them into the pangenome.

    This function processes system data from a single source, including both
    regular systems and their canonical representations if requested.

    Args:
        pangenome (Pangenome): Target the pangenome object to populate with systems.
        source_group (tables.Group): HDF5 group containing system tables for one source.
        models (Models): Model definitions associated with the systems.
        read_canonical (bool, optional): Whether to read canonical system representations.
            Defaults to True.
        disable_bar (bool, optional): Whether to disable the progress bar display.
            Defaults to False.

    Note:
        Systems are sorted by complexity before addition to ensure consistent ordering.
    """
    source = source_group._v_name
    logger.info(f"Reading systems from source: {source}")

    # Create mapping from functional unit names to model names for efficient lookup
    fu2model = {fu.name: model.name for model in models for fu in model.func_units}

    # Get table references
    system_table = source_group.systems
    unit_table = source_group.units

    def _read_system_unit(unit_row: tables.Table.row, units_dict: Dict[str, SystemUnit]) -> None:
        """
        Process a single unit row and add it to the units' dictionary.

        Args:
            unit_row (tables.Table.row): Row containing unit information.
            units_dict (Dict[str, SystemUnit]): Dictionary to store processed units.
        """
        unit_id = unit_row["ID"]

        # Create a new unit if not already processed
        if unit_id not in units_dict:
            fu_name = unit_row["name"].decode()
            model = models.get_model(fu2model[fu_name])
            unit = SystemUnit(functional_unit=model.get(fu_name), source=source)
            unit.ID = unit_id
            units_dict[unit_id] = unit
        else:
            unit = units_dict[unit_id]

        # Add gene family to the unit
        unit.add_family(
            pangenome.get_gene_family(unit_row["geneFam"].decode()),
            unit_row["metadata_source"].decode(),
            int(unit_row["metadata_id"])
        )

    def _read_system(
            sys_row: tables.Table.row,
            sys_dict: Dict[str, System],
            unit_dict: Dict[str, SystemUnit]
    ) -> System:
        """
        Process a single system row and add it to the systems' dictionary.

        Args:
            sys_row (tables.Table.row): Row containing system information.
            sys_dict (Dict[str, System]): Dictionary to store processed systems.
            unit_dict (Dict[str, SystemUnit]): Dictionary of available units.

        Returns:
            System: The processed system object.
        """
        sys_id = sys_row["ID"].decode()

        # Create a new system if not already processed
        if sys_id not in sys_dict:
            model = models.get_model(sys_row["name"].decode())
            sys = System(system_id=sys_id, model=model, source=source)
            sys_dict[sys_id] = sys
        else:
            sys = sys_dict[sys_id]

        # Add unit to the system
        sys.add_unit(unit_dict[sys_row["unit"]])
        return sys

    # Process regular systems
    systems = {}
    total_rows = system_table.nrows + unit_table.nrows

    with tqdm(total=total_rows, unit="line", desc="Reading pangenome systems", disable=disable_bar) as progress:
        # First pass: read all units
        units = {}
        for row in read_chunks(unit_table):
            _read_system_unit(row, units)
            progress.update()

        # Second pass: read systems and link to units
        for row in read_chunks(system_table):
            _read_system(row, systems, units)
            progress.update()

    logger.debug(f"Successfully read {len(systems)} systems from {source}")

    # Process canonical systems if requested
    if read_canonical:
        logger.debug(f"Reading canonical systems from {source}")

        canonic = {}
        canon_table = source_group.canonic
        canon_unit_table = source_group.canonic_units
        sys2canonical_table = source_group.system_to_canonical

        canonical_total = (canon_table.nrows + canon_unit_table.nrows +
                           sys2canonical_table.nrows)

        with tqdm(total=canonical_total, unit="line",
                  desc="Reading canonical systems", disable=disable_bar) as progress:

            # Read canonical units
            canon_units = {}
            for row in read_chunks(canon_unit_table):
                _read_system_unit(row, canon_units)
                progress.update()

            # Read canonical systems
            for row in read_chunks(canon_table):
                _read_system(row, canonic, canon_units)
                progress.update()

            # Link systems to their canonical representations
            for row in read_chunks(sys2canonical_table):
                system_id = row["system"].decode()
                canonical_id = row["canonic"].decode()
                systems[system_id].add_canonical(canonic[canonical_id])
                progress.update()

        logger.debug(f"Successfully read {len(canonic)} canonical systems from {source}")

    # Add systems to pangenome in sorted order for consistency
    # Sort by: canonical model complexity, system size (descending), family count (descending)
    logger.info(f"Integrating systems from {source} into pangenome...")

    sorted_systems = sorted(
        systems.values(),
        key=lambda x: (len(x.model.canonical), -len(x), -x.number_of_families)
    )

    # Add systems with progress tracking only in debug mode
    show_progress = not disable_bar and logger.level == logging.DEBUG
    for system in tqdm(sorted_systems, disable=not show_progress, desc="Adding systems to pangenome"):
        pangenome.add_system(system)

    system_count = pangenome.number_of_systems(source=source, with_canonical=False)
    logger.info(f"Successfully integrated {system_count} systems from {source}")


def read_systems(
        pangenome: Pangenome,
        h5f: tables.File,
        models: List[Models],
        sources: List[str],
        read_canonical: bool = False,
        disable_bar: bool = False
) -> Set[str]:
    """
    Read system information from all sources in the pangenome HDF5 file.

    Args:
        pangenome (Pangenome): Target pangenome object.
        h5f (tables.File): Open HDF5 file handle containing pangenome data.
        models (List[Models]): List of model definitions for each source.
        sources (List[str]): List of source identifiers to process.
        read_canonical (bool, optional): Whether to read canonical representations.
            Defaults to False.
        disable_bar (bool, optional): Whether to disable progress bars.
            Defaults to False.

    Returns:
        Set[str]: Combined set of metadata sources from all processed sources.

    Raises:
        ValueError: If the number of models doesn't match the number of sources.
    """
    if len(models) != len(sources):
        raise ValueError(f"Number of models ({len(models)}) must match number of sources ({len(sources)})")

    systems_group = h5f.root.systems
    metadata_sources = set()

    # Process each source sequentially
    for index, source in enumerate(sources):
        try:
            source_group = h5f.get_node(systems_group, source)

            # Accumulate metadata sources
            if hasattr(source_group._v_attrs, 'metadata_sources'):
                metadata_sources |= set(source_group._v_attrs.metadata_sources)

            # Read systems from this source
            read_systems_by_source(
                pangenome=pangenome,
                source_group=source_group,
                models=models[index],
                read_canonical=read_canonical,
                disable_bar=disable_bar
            )

            logger.debug(f"Completed processing source: {source}")

        except tables.NoSuchNodeError as e:
            logger.error(f"Source '{source}' not found in systems group: {e}")
            raise
        except Exception as e:
            logger.error(f"Error processing source '{source}': {e}")
            raise

    # Update pangenome status
    pangenome.status["systems"] = "Loaded"
    logger.info(f"Successfully loaded systems from {len(sources)} sources")

    return metadata_sources


def read_gene_families_info(
        pangenome: Pangenome,
        h5f: tables.File,
        information: bool = False,
        sequences: bool = False,
        disable_bar: bool = False
) -> None:
    """
    Read additional information about gene families from the HDF5 file.

    This function can read partition information and/or protein sequences
    for gene families already present in the pangenome.

    Args:
        pangenome (Pangenome): Target pangenome object containing gene families.
        h5f (tables.File): Open HDF5 file handle with gene family information.
        information (bool, optional): Whether to read partition information.
            Defaults to False.
        sequences (bool, optional): Whether to read protein sequences.
            Defaults to False.
        disable_bar (bool, optional): Whether to disable the progress bar.
            Defaults to False.

    Note:
        At least one of 'information' or 'sequences' should be True for this
        function to perform any meaningful work.
    """
    if not (information or sequences):
        logger.warning("Neither information nor sequences requested - no data will be read")
        return

    table = h5f.root.geneFamiliesInfo

    # Build a descriptive progress message
    components = []
    if information:
        components.append("information")
    if sequences:
        components.append("sequences")
    description = f"Reading gene families {' and '.join(components)}"

    # Process gene family information in chunks for memory efficiency
    chunk_size = 20000
    for row in tqdm(
            read_chunks(table, chunk=chunk_size),
            total=table.nrows,
            unit="gene family",
            desc=description,
            disable=disable_bar
    ):
        # Get or create a gene family
        family_name = row["name"].decode()
        fam = pangenome.get_gene_family(family_name)

        # Add requested information
        if information:
            fam.partition = row["partition"].decode()
        if sequences:
            fam.add_sequence(row["protein"].decode())

    # Update pangenome status based on what was loaded
    if information and h5f.root.status._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "Loaded"
        logger.debug("Gene family partition information loaded")

    if sequences and h5f.root.status._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "Loaded"
        logger.debug("Gene family sequences loaded")


def read_gene_families(
        pangenome: Pangenome,
        h5f: tables.File,
        disable_bar: bool = False
) -> None:
    """
    Read gene family associations from the HDF5 file.

    This function creates gene families and associates genes with them.
    If genome annotations are already loaded, it will link to existing
    gene objects; otherwise, it creates minimal gene objects.

    Args:
        pangenome (Pangenome): Target pangenome object.
        h5f (tables.File): Open HDF5 file handle containing gene family data.
        disable_bar (bool, optional): Whether to disable the progress bar.
            Defaults to False.
    """
    table = h5f.root.geneFamilies

    # Determine if we should link to existing annotations
    annotations_loaded = pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]
    logger.debug(f"Genome annotations loaded: {annotations_loaded}")

    # Process gene-family associations in chunks
    chunk_size = 20000
    for row in tqdm(
            read_chunks(table, chunk=chunk_size),
            total=table.nrows,
            unit="gene",
            desc="Associating genes with gene families",
            disable=disable_bar
    ):
        family_name = row["geneFam"].decode()
        gene_id = row["gene"].decode()

        # Get or create a gene family
        try:
            fam = pangenome.get_gene_family(name=family_name)
        except KeyError:
            # Create a new gene family with unique ID
            fam = GeneFamily(family_id=pangenome.max_fam_id, name=family_name)
            pangenome.add_gene_family(fam)

        # Get or create a gene object
        if annotations_loaded:
            # Link to existing annotated gene
            gene_obj = pangenome.get_gene(gene_id)
        else:
            # Create minimal gene object
            gene_obj = Gene(gene_id)

        # Associate gene with family
        fam.add(gene_obj)

    # Update pangenome status
    pangenome.status["genesClustered"] = "Loaded"
    logger.info("Successfully loaded gene family associations")


def read_spots(
        pangenome: Pangenome,
        h5f: tables.File,
        disable_bar: bool = False
) -> None:
    """
    Read genomic hotspots (spots) from the HDF5 file.

    Spots represent clusters of regions of genomic plasticity (RGPs)
    that occur in similar genomic contexts across multiple genomes.

    Args:
        pangenome (Pangenome): Target pangenome object.
        h5f (tables.File): Open HDF5 file handle with precomputed spots.
        disable_bar (bool, optional): Whether to disable the progress bar.
            Defaults to False.
    """
    table = h5f.root.spots
    spots = {}
    current_spot_id = None
    current_spot = None

    # Process spot data - regions are grouped by spot ID
    chunk_size = 20000
    for row in tqdm(
            read_chunks(table, chunk=chunk_size),
            total=table.nrows,
            unit="region",
            desc="Reading genomic spots",
            disable=disable_bar
    ):
        spot_id = int(row["spot"])

        # Create a new spot if we've moved to a different spot ID
        if current_spot_id != spot_id:
            current_spot_id = spot_id
            current_spot = spots.get(current_spot_id)

            if current_spot is None:
                current_spot = Spot(spot_id)
                spots[spot_id] = current_spot

        # Add a region to the current spot
        region_id = row["RGP"].decode()
        region = pangenome.get_region(region_id)
        current_spot.add(region)

    # Process spots and add to pangenome
    for spot in spots.values():
        # Generate family associations for the spot
        spot.spot_2_families()
        pangenome.add_spot(spot)

    # Update status
    pangenome.status["spots"] = "Loaded"
    logger.info(f"Successfully loaded {len(spots)} genomic spots")


def read_modules(
        pangenome: Pangenome,
        h5f: tables.File,
        disable_bar: bool = False
) -> None:
    """
    Read functional modules from the HDF5 file.

    Modules represent sets of gene families that are consistently
    found together and likely represent functional units.

    Args:
        pangenome (Pangenome): Target pangenome object.
        h5f (tables.File): Open HDF5 file handle with precomputed modules.
        disable_bar (bool, optional): Whether to disable the progress bar.
            Defaults to False.

    Raises:
        Exception: If gene families have not been loaded into the pangenome.
    """
    # Validate prerequisites
    if pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
        raise Exception(
            "Gene families must be loaded before reading modules. "
            "Please load gene families first."
        )

    table = h5f.root.modules
    modules = {}  # module_id -> Module object

    # Process module data in chunks
    chunk_size = 20000
    for row in tqdm(
            read_chunks(table, chunk=chunk_size),
            total=table.nrows,
            unit="association",
            desc="Reading functional modules",
            disable=disable_bar
    ):
        module_id = int(row["module"])
        family_name = row["geneFam"].decode()

        # Get or create module
        if module_id not in modules:
            modules[module_id] = Module(module_id)

        # Add gene family to module
        family = pangenome.get_gene_family(family_name)
        modules[module_id].add(family)

    # Add all modules to pangenome
    for module in modules.values():
        pangenome.add_module(module)

    # Update status
    pangenome.status["modules"] = "Loaded"
    logger.info(f"Successfully loaded {len(modules)} functional modules")


def read_pangenome(
        pangenome: Pangenome,
        annotation: bool = False,
        gene_families: bool = False,
        graph: bool = False,
        rgp: bool = False,
        spots: bool = False,
        gene_sequences: bool = False,
        modules: bool = False,
        metadata: bool = False,
        systems: bool = False,
        disable_bar: bool = False,
        **kwargs
) -> None:
    """
    Read a pangenome from its HDF5 file with specified components.

    This is the main function for loading pangenome data. It reads only
    the requested components and validates that they are available in the file.

    Args:
        pangenome (Pangenome): Target the pangenome object with the associated HDF5 file.
        annotation (bool, optional): Read genome annotations. Defaults to False.
        gene_families (bool, optional): Read gene family associations. Defaults to False.
        graph (bool, optional): Read gene neighborhood graph. Defaults to False.
        rgp (bool, optional): Read regions of genomic plasticity. Defaults to False.
        spots (bool, optional): Read genomic hotspots. Defaults to False.
        gene_sequences (bool, optional): Read gene DNA sequences. Defaults to False.
        modules (bool, optional): Read functional modules. Defaults to False.
        metadata (bool, optional): Read associated metadata. Defaults to False.
        systems (bool, optional): Read biological systems. Defaults to False.
        disable_bar (bool, optional): Disable all progress bars. Defaults to False.
        **kwargs: Additional parameters

    Raises:
        FileNotFoundError: If pangenome has no associated HDF5 file.
        ValueError: If requested data is not available in the file.
        AttributeError: If required graph/spots/modules data is missing.
        KeyError: If required metadata is not present.
    """
    # Validate input
    if pangenome.file is None:
        raise FileNotFoundError(
            "The pangenome object must have an associated HDF5 file. "
            "Use pangenome.add_file(path) to set the file path."
        )

    logger.info(f"Reading pangenome from: {pangenome.file}")

    # Open HDF5 file for reading
    with tables.open_file(pangenome.file, "r") as h5f:
        status_attrs = h5f.root.status._v_attrs

        # Read annotations if requested
        if annotation:
            if status_attrs.genomesAnnotated:
                logger.info("Reading pangenome annotations...")
                read_annotation(pangenome, h5f, disable_bar=disable_bar)
            else:
                raise ValueError(
                    f"Pangenome in '{pangenome.file}' has not been annotated "
                    "or has been improperly filled"
                )

        # Read gene sequences if requested
        if gene_sequences:
            if status_attrs.geneSequences:
                logger.info("Reading pangenome gene DNA sequences...")
                read_gene_sequences(pangenome, h5f, disable_bar=disable_bar)
            else:
                raise ValueError(
                    f"Pangenome in '{pangenome.file}' does not contain gene sequences "
                    "or has been improperly filled"
                )

        # Read gene families if requested
        if gene_families:
            if status_attrs.genesClustered:
                logger.info("Reading pangenome gene families...")
                read_gene_families(pangenome, h5f, disable_bar=disable_bar)

                # Read additional gene family information if requested
                if kwargs.get("gene_families_info") or kwargs.get("gene_families_sequences"):
                    info_components = []
                    if kwargs.get("gene_families_info"):
                        info_components.append("partition info")
                    if kwargs.get("gene_families_sequences"):
                        info_components.append("sequences")

                    logger.info(f"Reading gene family {' and '.join(info_components)}...")
                    read_gene_families_info(
                        pangenome, h5f,
                        information=kwargs.get("gene_families_info", False),
                        sequences=kwargs.get("gene_families_sequences", False),
                        disable_bar=disable_bar
                    )
            else:
                raise ValueError(
                    f"Pangenome in '{pangenome.file}' does not contain gene families "
                    "or has been improperly filled"
                )

        # Read gene neighborhood graph if requested
        if graph:
            if status_attrs.NeighborsGraph:
                logger.info("Reading gene neighborhood graph...")
                read_graph(pangenome, h5f, disable_bar=disable_bar)
            else:
                raise AttributeError(
                    f"Pangenome in '{pangenome.file}' does not contain graph information "
                    "or has been improperly filled"
                )

        # Read regions of genomic plasticity if requested
        if rgp:
            if status_attrs.predictedRGP:
                logger.info("Reading regions of genomic plasticity...")
                read_rgp(pangenome, h5f, disable_bar=disable_bar)
            else:
                raise AttributeError(
                    f"Pangenome in '{pangenome.file}' does not contain RGP information "
                    "or has been improperly filled"
                )

        # Read genomic spots if requested
        if spots:
            if status_attrs.spots:
                logger.info("Reading genomic spots...")
                start_time = time.time()
                read_spots(pangenome, h5f, disable_bar=disable_bar)
                logger.debug(f"Loading spots took: {time.time() - start_time:.2f} seconds")
            else:
                raise AttributeError(
                    f"Pangenome in '{pangenome.file}' does not contain spots information "
                    "or has been improperly filled"
                )

        # Read functional modules if requested
        if modules:
            if status_attrs.modules:
                logger.info("Reading functional modules...")
                read_modules(pangenome, h5f, disable_bar=disable_bar)
            else:
                raise AttributeError(
                    f"Pangenome in '{pangenome.file}' does not contain modules information "
                    "or has been improperly filled"
                )

        # Read biological systems if requested
        if systems:
            # Systems reading requires metadata, so enable it
            metadata = True

            logger.info("Reading biological systems...")
            metadata_sources = read_systems(
                pangenome=pangenome,
                h5f=h5f,
                models=kwargs["models"],
                sources=kwargs["systems_sources"],
                read_canonical=kwargs.get("read_canonical", False),
                disable_bar=disable_bar
            )

            # Update metadata parameters for systems
            if "meta_sources" in kwargs:
                kwargs["metatypes"].add("families")
                kwargs["meta_sources"] |= metadata_sources
            else:
                kwargs["metatypes"] = {"families"}
                kwargs["meta_sources"] = metadata_sources

        # Read metadata if requested
        if metadata:
            metatypes = kwargs.get("metatypes", set())
            meta_sources = kwargs.get("meta_sources", set())
            sources = kwargs.get("sources", set())

            for metatype in metatypes:
                if status_attrs.metadata:
                    metastatus = h5f.root.status._f_get_child("metastatus")
                    metasources = h5f.root.status._f_get_child("metasources")

                    # Determine which sources to read for this metatype
                    available_sources = set(metasources._v_attrs[metatype])
                    requested_sources = available_sources & sources

                    if meta_sources:
                        requested_sources &= meta_sources

                    # Read metadata if available and requested
                    if metastatus._v_attrs[metatype] and requested_sources:
                        logger.info(
                            f"Reading {metatype} metadata from sources: {requested_sources}"
                        )
                        read_metadata(
                            pangenome, h5f, metatype, requested_sources,
                            disable_bar=disable_bar
                        )
                else:
                    raise KeyError(
                        f"Pangenome in '{pangenome.file}' does not contain "
                        f"metadata for {metatype}"
                    )

    logger.info("Pangenome reading completed successfully")


def check_pangenome_info(
        pangenome: Pangenome,
        need_families_info: bool = False,
        need_families_sequences: bool = False,
        need_systems: bool = False,
        models: Optional[List[Models]] = None,
        systems_sources: Optional[List[str]] = None,
        read_canonical: bool = False,
        disable_bar: bool = False,
        **kwargs
) -> None:
    """
    Determine and load required pangenome information based on analysis needs.

    This function analyzes what information is needed and automatically loads
    the required components from the pangenome file.

    Args:
        pangenome (Pangenome): Target pangenome object.
        need_families_info (bool, optional): Whether gene family partition info is needed.
            Defaults to False.
        need_families_sequences (bool, optional): Whether gene family sequences are needed.
            Defaults to False.
        need_systems (bool, optional): Whether biological systems are needed.
            Defaults to False.
        models (Optional[List[Models]], optional): Model definitions for systems.
            Required if need_systems=True.
        systems_sources (Optional[List[str]], optional): System source identifiers.
            Required if need_systems=True.
        read_canonical (bool, optional): Whether to read canonical system representations.
            Defaults to False.
        disable_bar (bool, optional): Whether to disable progress bars.
            Defaults to False.
        **kwargs: Additional parameters
    Raises:
        AssertionError: If systems are requested but required parameters are missing.
        ValueError: If gene families info/sequences are requested without gene families.
    """
    # Get base information requirements
    need_info = get_need_info(pangenome, **kwargs)
    logger.debug(f"Base information requirements: {need_info}")

    # Validate gene family information requirements
    if need_families_info or need_families_sequences:
        if not need_info.get('gene_families', False):
            raise ValueError(
                "Gene families must be loaded to access family information or sequences. "
                "Set gene_families=True in your request."
            )

    # Add gene family specific requirements
    need_info["gene_families_info"] = need_families_info
    need_info["gene_families_sequences"] = need_families_sequences

    # Handle systems requirements
    if need_systems:
        if models is None or systems_sources is None:
            raise AssertionError(
                "Both 'models' and 'systems_sources' parameters are required "
                "when need_systems=True"
            )

        need_info.update({
            "systems": True,
            "models": models,
            "systems_sources": systems_sources,
            "read_canonical": read_canonical
        })

    logger.debug(f"Information requirements determined: {need_info}")

    # Load required information if any is needed
    if any(need_info.values()):
        logger.info("Loading required pangenome information...")
        read_pangenome(pangenome, disable_bar=disable_bar, **need_info)
    else:
        logger.info("No additional information loading required")


def load_pangenome(
        name: str,
        path: Path,
        taxid: int,
        need_info: Dict[str, bool],
        check_function: Optional[Callable[[Pangenome, ...], None]] = None,
        disable_bar: bool = False,
        **kwargs
) -> Pangenome:
    """
    Load a single pangenome from the file with specified information requirements.

    This function creates a new pangenome object, associates it with the given
    HDF5 file, performs optional validation, and loads the requested information.

    Args:
        name (str): Descriptive name for the pangenome.
        path (Path): Path to the pangenome HDF5 file.
        taxid (int): NCBI taxonomic identifier for the pangenome.
        need_info (Dict[str, bool]): Dictionary specifying what information to load.
            Keys can include: 'annotation', 'gene_families', 'graph', 'rgp',
            'spots', 'gene_sequences', 'modules', 'metadata', 'systems', etc.
        check_function (Optional[Callable], optional): Custom validation function
            to run before loading information. Should raise exceptions on failure.
        disable_bar (bool, optional): Whether to disable progress bars.
            Defaults to False.
        **kwargs: Additional parameters passed to check_function and
            check_pangenome_info.

    Returns:
        Pangenome: Fully loaded pangenome object with requested information.

    Raises:
        FileNotFoundError: If the specified path does not exist.
        Exception: If check_function validation fails or loading encounters errors.
    """
    start_time = time.time()
    logger.info(f"Starting to load pangenome: {name}")

    # Validate file exists
    if not path.exists():
        raise FileNotFoundError(f"Pangenome file not found: {path}")

    # Create and configure the pangenome object
    pangenome = Pangenome(name=name, taxid=taxid)
    pangenome.add_file(path)

    # Run optional validation
    if check_function is not None:
        logger.debug(f"Running validation function for pangenome: {name}")
        try:
            check_function(pangenome, **kwargs)
        except Exception as error:
            logger.error(f"Validation failed for pangenome {name}: {error}")
            raise error

    # Load requested information
    try:
        check_pangenome_info(pangenome, disable_bar=disable_bar, **need_info)
    except Exception as error:
        logger.error(f"Failed to load information for pangenome {name}: {error}")
        raise error

    # Log completion
    elapsed_time = time.time() - start_time
    logger.info(f"Pangenome '{name}' loaded successfully in {elapsed_time:.2f} seconds")

    return pangenome


def load_pangenomes(
        pangenome_list: Path,
        need_info: Dict[str, bool],
        check_function: Optional[Callable] = None,
        max_workers: int = 1,
        lock: Optional[Lock] = None,
        disable_bar: bool = False,
        **kwargs
) -> Pangenomes:
    """
    Load multiple pangenomes in parallel from a configuration file.

    This function provides efficient parallel loading of multiple pangenomes
    using a thread pool. It reads pangenome specifications from a TSV file
    and loads each pangenome with the same information requirements.

    Args:
        pangenome_list (Path): Path to TSV file containing pangenome specifications.
            Expected format: name\tpath\ttaxid (tab-separated values).
        need_info (Dict[str, bool]): Information requirements applied to all pangenomes.
            See load_pangenome documentation for available keys.
        check_function (Optional[Callable], optional): Validation function applied
            to each pangenome before loading. Defaults to None.
        max_workers (int, optional): Maximum number of concurrent loading threads.
            Defaults to 1 (sequential loading).
        lock (Optional[Lock], optional): Multiprocessing lock for thread synchronization.
            If None, a new lock will be created. Defaults to None.
        disable_bar (bool, optional): Whether to disable the main progress bar.
            Individual pangenome progress bars are always disabled in parallel mode.
            Defaults to False.
        **kwargs: Additional parameters passed to each pangenome's load_pangenome call.

    Returns:
        Pangenomes: Collection containing all successfully loaded pangenomes.

    Raises:
        FileNotFoundError: If the pangenome_list file does not exist.
        ValueError: If the pangenome_list file format is invalid.
        Exception: If any pangenome fails to load (will stop all loading).

    Note:
        - Progress bars for individual pangenomes are disabled during parallel loading
          to avoid output conflicts.
        - All pangenomes must load successfully; if any fails, the entire operation fails.
        - The lock parameter is used to synchronize access to shared resources during
          parallel execution.
    """
    start_time = time.time()
    logger.info(f"Starting batch pangenome loading from: {pangenome_list}")

    # Validate and parse input file
    try:
        pan_to_path = check_tsv_sanity(pangenome_list)
    except Exception as error:
        logger.error(f"Failed to parse pangenome list file: {error}")
        raise

    if not pan_to_path:
        logger.warning("No pangenomes found in the input file")
        return Pangenomes()

    logger.info(f"Found {len(pan_to_path)} pangenomes to load")

    # Initialize result container and synchronization
    pangenomes = Pangenomes()
    if lock is None:
        lock = Lock()

    # Configure parallel execution
    # Disable individual progress bars during parallel loading to avoid conflicts
    individual_disable_bar = max_workers > 1 or disable_bar

    with ThreadPoolExecutor(
            max_workers=max_workers,
            initializer=init_lock,
            initargs=(lock,)
    ) as executor:

        # Submit all loading tasks
        with tqdm(
                total=len(pan_to_path),
                unit='pangenome',
                desc="Loading pangenomes",
                disable=disable_bar
        ) as progress:

            futures = []
            for pangenome_name, pangenome_info in pan_to_path.items():
                # Extract pangenome information
                pangenome_path = pangenome_info["path"]
                pangenome_taxid = pangenome_info["taxid"]

                # Submit loading task
                future = executor.submit(
                    load_pangenome,
                    name=pangenome_name,
                    path=pangenome_path,
                    taxid=pangenome_taxid,
                    need_info=need_info,
                    check_function=check_function,
                    disable_bar=individual_disable_bar,
                    **kwargs
                )

                # Set up a progress callback
                future.add_done_callback(lambda f: progress.update())
                futures.append((future, pangenome_name))

            # Collect results as they complete
            for future, pangenome_name in futures:
                try:
                    # Get the result (this will raise any exceptions that occurred)
                    pangenome_result = future.result()

                    # Thread-safe addition to results
                    with lock:
                        pangenomes.add(pangenome_result)

                except Exception as error:
                    logger.error(f"Failed to load pangenome '{pangenome_name}': {error}")
                    # Cancel remaining futures and re-raise
                    for remaining_future, _ in futures:
                        remaining_future.cancel()
                    raise error

    # Log completion statistics
    elapsed_time = time.time() - start_time
    logger.info(
        f"Successfully loaded {len(pangenomes)} pangenomes in "
        f"{elapsed_time:.2f} seconds"
    )

    if max_workers > 1:
        avg_time_per_pangenome = elapsed_time / len(pangenomes) if pangenomes else 0
        logger.debug(
            f"Average loading time per pangenome: {avg_time_per_pangenome:.2f} seconds "
            f"(with {max_workers} workers)"
        )

    return pangenomes
