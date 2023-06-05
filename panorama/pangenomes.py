#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Generator, Iterable, List, Set, Union

# install libraries
from ppanggolin.pangenome import Pangenome as Pan

# local libraries
from panorama.system import System
from panorama.geneFamily import GeneFamily
from panorama.region import Module


class Pangenome(Pan):
    """
    This is a class representing pangenome based on PPanGGOLLiN class. It is used as a basic unit for all the analysis
    to access to the different elements of your pangenome, such as organisms, contigs, genes or gene families.
    This class provide some more methods needed to analyse pangenome.

    :param name: Name of the pangenome
    """

    def __init__(self, name, taxid: int = None):
        """Constructor method.
        """
        super().__init__()
        self._system_getter = {}
        self._name2system = {}
        self._max_id_system = 0
        self.name = name
        self.taxid = taxid
        self.status.update({"systems": 'No',
                            "systems_sources": set()})

    def __str__(self):
        return self.name

    def add_file(self, pangenome_file: str):
        """Links an HDF5 file to the pan. If needed elements will be loaded from this file,
        and anything that is computed will be saved to this file when
        :func:`ppanggolin.formats.writeBinaries.writePangenome` is called.
        :param pangenome_file: A string representing the filepath to the hdf5 pan file
        to be either used or created
        """
        from panorama.format.read_binaries import get_status
        # importing on call instead of importing on top to avoid cross-reference problems.
        get_status(self, pangenome_file)
        self.file = pangenome_file

    @property
    def gene_families(self) -> List[GeneFamily]:
        """returns all the gene families in the pangenome

        :return: list of gene families
        """
        return super().gene_families

    def _create_gene_family(self, name: str) -> GeneFamily:
        """Creates a gene family object with the given `name`

        :param name: the name to give to the gene family. Must not exist already.
        :return: the created GeneFamily object
        """

        new_fam = GeneFamily(family_id=self.max_fam_id, name=name)
        self.max_fam_id += 1
        self._famGetter[new_fam.name] = new_fam
        return new_fam

    def get_gene_family(self, name: str) -> Union[GeneFamily, None]:
        """ Get the gene family by his name in the pangenome

        :param name: name of the gene family to get

        :return: The desired gene family

        :raise KeyError: Gene family doesn't exist in pangenome
        :raise Exception: Manage unexpected error
        """
        try:
            fam = super().get_gene_family(name)
        except KeyError:
            raise KeyError(f"{name} is not find in pangenome gene families")
        except Exception:
            raise Exception(f"Unexpected problems to get {name} gene familiy in pangenome")
        else:
            return fam

    @property
    def systems(self) -> Generator[System, None, None]:
        """Get all systems in pangenome
        """
        for system in self._system_getter.values():
            yield system

    @property
    def systems_sources(self) -> Set[str]:
        sources = set()
        for system in self.systems:
            sources.add(system.source)
        return sources

    def get_system(self, system_id: str) -> System:
        try:
            system = self._system_getter[system_id]
        except KeyError:
            uncanonical_id = system_id.split('.')[0]
            find_in_canonical = False
            if len(uncanonical_id) > 0:
                uncanonical_system = self.get_system(uncanonical_id)
                for canonical in uncanonical_system.canonical:
                    if canonical.ID == system_id:
                        find_in_canonical = True
                        return canonical
                if not find_in_canonical:
                    raise KeyError(f"There is no system with ID = {system_id} in pangenome")
            else:
                raise KeyError(f"There is no system with ID = {system_id} in pangenome")
        else:
            return system

    def get_system_by_source(self, source: str) -> Generator[System, None, None]:
        for system in self.systems:
            if system.source == source:
                yield system

    def add_system(self, system: System):
        """Add a detected system in the pangenome

        :param system: Detected system that will be added
        """
        same_sys = False
        canonical_systems = []
        drop_sys_key = []
        for system_in in self.get_system_by_source(system.source):
            if system_in.name == system.name and system_in.gene_families.issubset(system.gene_families):
                # A system with this name already exist and system in pangenome is subset of new system
                system.ID = system_in.ID
                self._system_getter[system.ID] = system
                same_sys = True
            elif system.name in system_in.canonical_models():
                # System in pangenome is a canonical system for new system
                if len(system.gene_families.intersection(system_in.gene_families)) > 0:
                    canonical_systems.append(system_in)
                    drop_sys_key.append(system_in.ID)
            elif system_in.name in system.canonical_models():
                # New system is a canonical system for a system in pangenome
                if len(system.gene_families.intersection(system_in.gene_families)) > 0:
                    system_in.add_canonical(system)
                    same_sys = True

        if not same_sys:
            self._max_id_system += 1
            system.ID = str(self._max_id_system)
            self._system_getter[system.ID] = system
            for canonical_system in canonical_systems:
                system.add_canonical(canonical_system)
        self._system_getter = {sys_id: sys for sys_id, sys in self._system_getter.items() if sys_id not in drop_sys_key}

    def number_of_systems(self, source: str = None, with_canonical: bool = True) -> int:
        """Get the number of systems in the pangenomes"""
        nb_systems = 0
        systems = self.systems if source is None else self.get_system_by_source(source)
        for system in systems:
            nb_systems += 1 + len(system.canonical) if with_canonical else 1
        return nb_systems
    
    def add_modules(self, modules: Iterable[Module]):
        super().add_modules({Module(module_id=module.ID, families=module.families) for module in modules})


class Pangenomes:
    """
    This class represente a group of pangenome object.
    """

    def __init__(self):
        """Constructor method
        """
        self._pangenomes_getter = dict()

    def __len__(self):
        return len(self._pangenomes_getter)

    def __iter__(self) -> Generator[Pangenome, None, None]:
        for pangenome in self._pangenomes_getter.values():
            yield pangenome

    def add_pangenome(self, pangenome: Pangenome):
        """ Add a pangenome object

        :param pangenome: Pangenome object
        :return:
        """
        self._pangenomes_getter[pangenome.name] = pangenome

    def add_list_pangenomes(self, pangenomes_list: List[Pangenome]):
        for pangenome in pangenomes_list:
            self.add_pangenome(pangenome)


    def to_list(self) -> List[Pangenome]:
        return list(self.__iter__())

    def to_set(self) -> Set[Pangenome]:
        return set(self.__iter__())

    def get_pangenome(self, name: str):
        return self._pangenomes_getter[name]
