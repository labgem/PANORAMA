#!/usr/bin/env python3
# coding: utf8

# default libraries
import logging
from pathlib import Path
from typing import Dict, Generator, List, Set, Union
from collections import defaultdict

# install libraries
import pandas as pd
from tqdm import tqdm
from ppanggolin.pangenome import Pangenome as Pan
from ppanggolin.pangenome import GeneFamily as Fam

# local libraries
from panorama.systems.system import System
from panorama.geneFamily import GeneFamily, Akin
from panorama.region import Spot, ConservedSpots


class Pangenome(Pan):
    """This is a class representing pangenome based on PPanGGOLLiN class.
    It is used as a basic unit for all the analysis to access to the different elements
    of your pangenome, such as organisms, contigs, genes or gene families.
    This class provides some more methods needed to analyze pangenome.

    Args:
        name: Name of the pangenome
    """

    def __init__(self, name, taxid: int = None):
        """Constructor method.
        """
        super().__init__()
        self._system_getter = {}
        self._name2system = {}
        self._max_id_system = 0
        self._systems_sources = None
        self._systems_sources2metada_sources = None
        self.name = name
        self.taxid = taxid
        self.status.update({"systems": 'No', "systems_sources": set()})

    def __str__(self):
        return self.name

    def add_file(self, pangenome_file: Path, check_version: bool = True):
        """Links an HDF5 file to the pan.

        If needed elements will be loaded from this file,
        and anything that is computed will be saved to this file when
        :func:`ppanggolin.formats.writeBinaries.writePangenome` is called.

        Args:
            pangenome_file: A string representing the filepath to the hdf5 pan file
                to be either used or created
            check_version: Check ppanggolin version of the pangenome file to be compatible with the current version of ppanggolin being used.
        Raises:
            AssertionError: If the `pangenome_file` is not an instance of the Path class
            TypeError: If the `pangenome_file` is not a HDF5 format file
        """
        assert isinstance(pangenome_file, Path), "pangenome file should be a Path object type"
        from tables import is_hdf5_file
        from ppanggolin.utils import check_version_compatibility
        from panorama.format.read_binaries import get_status
        # importing on call instead of importing on top to avoid cross-reference problems.
        if not is_hdf5_file(pangenome_file):
            raise TypeError("Pangenome file should be an HDF5 file type")
        get_status(self, pangenome_file)

        check_version_compatibility(self.status["ppanggolin_version"])

        self.file = pangenome_file.absolute().as_posix()

    @property
    def gene_families(self) -> Generator[GeneFamily, None, None]:
        """Returns all the gene families in the pangenome

        Returns:
            Generator[GeneFamily, None, None]: Generator of gene families
        """
        return super().gene_families

    def _create_gene_family(self, name: str) -> GeneFamily:
        """Creates a gene family object with the given `name`

        Args:
            name: The name to give to the gene family. Must not exist already.

        Returns:
            GeneFamily: The created GeneFamily object
        """
        new_fam = GeneFamily(family_id=self.max_fam_id, name=name)
        self.max_fam_id += 1
        self._fam_getter[new_fam.name] = new_fam
        return new_fam

    def add_gene_family(self, family: Union[GeneFamily, Fam]):
        """Adds a gene family to the pangenome

        Args:
            family (Union[GeneFamily, Fam]): GeneFamily object to add
        """
        assert isinstance(family, (GeneFamily, Fam)), "Family must be a GeneFamily from PANORAMA or PPanGGOLiN"
        if not isinstance(family, GeneFamily):
            family = GeneFamily.recast(family)
        super().add_gene_family(family)

    def get_gene_family(self, name: str) -> Union[GeneFamily, None]:
        """Get the gene family by its name in the pangenome

        Args:
            name: Name of the gene family to get

        Returns:
            Union[GeneFamily, None]: The desired gene family
        """
        return super().get_gene_family(name)

    def add_spot(self, spot: Spot):
        super().add_spot(spot)
        spot.pangenome = self

    @property
    def systems(self) -> Generator[System, None, None]:
        """Get all systems in the pangenome

        Yields:
            Generator[System, None, None]: Generator of systems
        """
        for system in self._system_getter.values():
            yield system

    def _get_systems_sources_and_related_metadata_sources(self):
        sources = set()
        sys_sources2meta_sources = defaultdict(set)
        for system in self.systems:
            sources.add(system.source)
            sys_sources2meta_sources[system.source] |= system.annotation_sources()
        self._systems_sources = sources
        self._systems_sources2metada_sources = sys_sources2meta_sources

    @property
    def systems_sources(self) -> Set[str]:
        """Get sources of all systems in the pangenome

        Returns:
            Set[str]: Set of system sources
        """
        if self._systems_sources is not None:
            return self._systems_sources
        else:
            self._get_systems_sources_and_related_metadata_sources()
            return self.systems_sources

    def systems_sources_to_metadata_source(self) -> Dict[str, Set[str]]:
        """Get metadata sources related to system sources

        Returns:
            Dict[str, Set[str]]: System source as key linked to their metadata sources as value
        """
        if self._systems_sources2metada_sources is not None:
            return self._systems_sources2metada_sources
        else:
            self._get_systems_sources_and_related_metadata_sources()
            return self.systems_sources_to_metadata_source()

    def get_system(self, system_id: str) -> System:
        """Get a system by its ID in the pangenome

        Args:
            system_id: ID of the system to get

        Returns:
            System: The desired system

        Raises:
            KeyError: If the system doesn't exist in the pangenome
        """
        try:
            system = self._system_getter[system_id]
        except KeyError:
            uncanonical_id = system_id.split('.')[0]
            find_in_canonical = False
            if len(uncanonical_id) > 0:
                uncanonical_system = self.get_system(uncanonical_id)
                for canonical in uncanonical_system.canonical:
                    if canonical.ID == system_id:
                        return canonical
                if not find_in_canonical:
                    raise KeyError(f"There is no system with ID = {system_id} in pangenome")
                else:
                    # You should not arrive here because if `find_in_canonical` is True the function return a value
                    raise Exception("Something unexpected happened here. Please report an issue on our GitHub")
            else:
                raise KeyError(f"There is no system with ID = {system_id} in pangenome")
        else:
            return system

    def get_system_by_source(self, source: str) -> Generator[System, None, None]:
        """Retrieve systems by their source.

        Args:
            source (str): Source identifier.

        Yields:
            Generator[System, None, None]: Systems with the given source.
        """
        for system in self.systems:
            if system.source == source:
                yield system

    def add_system(self, system: System):
        """Add a detected system to the pangenome.

        Args:
            system (System): Detected system to be added.

        TODO:
            Merge systems.
        """
        same_sys = False
        canonical_systems = []
        drop_sys_key = []
        for system_in in self.get_system_by_source(system.source):
            if system_in.name == system.name:
                if system_in.is_subset(system):
                    # A system with this name already exist and system in pangenome is subset of new system
                    system.ID = system_in.ID
                    self._system_getter[system.ID] = system
                    same_sys = True
                    for family in system.difference(system_in):
                        family.add_system(system)
                elif system_in.is_superset(system):
                    same_sys = True
            elif system.name in system_in.canonical_models():
                # System in pangenome is a canonical system for new system
                if len(system.intersection(system_in)) > 0:
                    canonical_systems.append(system_in)
                    drop_sys_key.append(system_in.ID)
            elif system_in.name in system.canonical_models():
                # New system is a canonical system for a system in pangenome
                if len(system.intersection(system_in)) > 0:
                    system_in.add_canonical(system)
                    self._max_id_system += 1
                    system.ID = str(self._max_id_system)
                    same_sys = True
                    for family in system.families:
                        family.add_system(system)

        if not same_sys:
            self._max_id_system += 1
            system.ID = str(self._max_id_system)
            self._system_getter[system.ID] = system
            for canonical_system in canonical_systems:
                system.add_canonical(canonical_system)
            for family in system.families:
                family.add_system(system)
        self._system_getter = {sys_id: sys for sys_id, sys in self._system_getter.items() if sys_id not in drop_sys_key}

    def number_of_systems(self, source: str = None, with_canonical: bool = True) -> int:
        """Get the number of systems in the pangenome.

        Args:
            source (str, optional): Source identifier. Defaults to None.
            with_canonical (bool, optional): Include canonical systems. Defaults to True.

        Returns:
            int: Number of systems.
        """
        nb_systems = 0
        systems = self.systems if source is None else self.get_system_by_source(source)
        for system in systems:
            nb_systems += 1 + len(system.canonical) if with_canonical else 1
        return nb_systems


class Pangenomes:
    """A collection of pangenome objects."""

    def __init__(self):
        """Initialize an empty collection of pangenomes."""
        self._pangenomes_getter = {}
        self._clusters = {}
        self._families2pangenome = {}
        self._conserved_spots_getter = {}

    def __len__(self):
        """Get the number of pangenomes in the collection."""
        return len(self._pangenomes_getter)

    def __iter__(self) -> Generator[Pangenome, None, None]:
        """Iterate over the pangenomes in the collection."""
        for pangenome in self._pangenomes_getter.values():
            yield pangenome

    def add(self, pangenome: Pangenome):
        """Add a pangenome object to the collection.

        Args:
            pangenome (Pangenome): The pangenome object to add.
        """
        self._pangenomes_getter[pangenome.name] = pangenome

    def to_list(self) -> List[Pangenome]:
        """Convert the collection to a list of pangenomes.

        Returns:
            List[Pangenome]: A list of pangenome objects.
        """
        return list(self.__iter__())

    def to_set(self) -> Set[Pangenome]:
        """Convert the collection to a set of pangenomes.

        Returns:
            Set[Pangenome]: A set of pangenome objects.
        """
        return set(self.__iter__())

    def get(self, name: str):
        """Retrieve a pangenome object from the collection by its name.

        Args:
            name (str): The name of the pangenome to retrieve.

        Returns:
            Pangenome: The pangenome object with the specified name.
        """
        return self._pangenomes_getter[name]

    def mk_families_to_pangenome(self, check_duplicate_names: bool = True):
        """
        Fill the families2pangenome dictionary to know from which pangenome the family belongs to

        Args:
            check_duplicate_names: Flag to return an error if families name is duplicated between pangenome.

        Raises:
            KeyError: If there is a duplicate family names
        """
        for pangenome in self:
            for family in pangenome.gene_families:
                try:
                    duplicate_pangenome = self._families2pangenome[family.name]
                except KeyError:
                    self._families2pangenome[family.name] = pangenome.name
                else:
                    if check_duplicate_names:
                        logging.getLogger("PANORAMA").error(f"Duplicate family name is {family.name} between "
                                                            f"{pangenome.name} and {duplicate_pangenome}")
                        raise KeyError("There is duplicate names of gene families between your pangenomes."
                                       "In your command it could affect the results. Look at the documentation or "
                                       "post an issue in our GitHub to manage this situation.")
                    else:
                        logging.getLogger("PANORAMA").debug("Duplicate family names between pangenomes")
                        # "Be really careful with what you do"
                        duplicate_family = self.get(duplicate_pangenome).get_gene_family(family.name)
                        duplicate_family.name = f"{duplicate_pangenome}_{family.name}"
                        family.name = f"{pangenome.name}_{family.name}"

    def get_family(self, name: str, check_duplicate_names: bool = True) -> GeneFamily:
        """
        Get a family in the pangenomes

        Args:
            name: name of the gene family
            check_duplicate_names: Flag to raise an error if duplicate family names between pangenomes.

        Returns:
            The gene family with the given name
        """
        if len(self._families2pangenome) == 0:
            self.mk_families_to_pangenome(check_duplicate_names)
        return self.get(self._families2pangenome[name]).get_gene_family(name)

    def add_cluster(self, cluster: Akin) -> None:
        """
        Add a cluster of similar gene families between pangenomes

        Args:
            cluster: A set of akin gene families

        Raises:
            KeyError: If there is already an Akin object with this ID
        """
        try:
            _ = self._clusters[cluster.ID]
        except KeyError:
            self._clusters[cluster.ID] = cluster
        else:
            raise KeyError(f"Cluster with ID {cluster.ID} already exist in pangenomes")

    def read_clustering(self, clustering: Union[Path, pd.DataFrame], disable_bar: bool = False):
        """
        Read clustering result from panorama

        Args:
            clustering: Clustering result
            disable_bar: Flag to disable progress bar (default: False)
        """
        from panorama.alignment.cluster import clust_col_names

        logging.getLogger("PANORAMA").info("Reading clustering...")
        if isinstance(clustering, Path):
            logging.getLogger("PANORAMA").debug(f"Reading clustering from {clustering}")
            clustering = pd.read_csv(clustering, sep="\t", names=clust_col_names, header=0)
        if clustering.shape[1] != len(clust_col_names) or clustering.columns.tolist() != clust_col_names:
            raise pd.errors.ParserError("given Clustering has inconsistent number of columns or names")
        clustering = clustering.sort_values(by=clust_col_names[0])
        lines = clustering.iterrows()
        line = next(lines)[1]
        cluster_id = line[clust_col_names[0]]
        referent = self.get_family(line[clust_col_names[1]])
        gene_families = set()
        stop = False
        with tqdm(total=clustering.shape[0], unit="line", disable=disable_bar) as pbar:
            while not stop:
                try:
                    line = next(lines)[1]
                except StopIteration:
                    stop = True
                else:
                    curr_cluster_id = line[clust_col_names[0]]
                    if curr_cluster_id != cluster_id:
                        self.add_cluster(Akin(cluster_id, referent, *gene_families))
                        cluster_id = curr_cluster_id
                        gene_families = set()
                    ref_name = line[clust_col_names[1]]
                    in_clust_name = line[clust_col_names[2]]
                    if ref_name == in_clust_name:
                        referent = self.get_family(ref_name)
                    else:
                        gene_families.add(self.get_family(in_clust_name))
                finally:
                    pbar.update()
        self.add_cluster(Akin(cluster_id, referent, *gene_families))

    @property
    def conserved_spots(self) -> Generator[ConservedSpots, None, None]:
        """Generator of conserved spots between pangenomes
        Yields:
            ConservedSpots: a set of spots conserved between pangenomes
        """
        for conserved_spot in self._conserved_spots_getter.values():
            yield conserved_spot

    def add_conserved_spots(self, conserved_spots: ConservedSpots):
        """
        Add a set of conserved spots between pangenomes

        Args:
            conserved_spots: Conserved spots object

        Raises:
            KeyError: if conserved_spots identifier already exist in pangenomes.
        """
        try:
            _ = self._conserved_spots_getter[conserved_spots.ID]
        except KeyError:
            self._conserved_spots_getter[conserved_spots.ID] = conserved_spots
        else:
            raise KeyError(f"Conserved spots {conserved_spots.ID} already exists between pangenomes")

    @property
    def number_of_conserved_spots(self) -> int:
        """Get the number of conserved spots

        Returns:
            Number of conserved spots
        """
        return len(self._conserved_spots_getter)
