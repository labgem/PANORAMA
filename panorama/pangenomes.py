#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Generator, List, Set, Union

# install libraries
from ppanggolin.pangenome import Pangenome as Pan

# local libraries
from annotation import Annotation
from system import System
from panorama.geneFamily import GeneFamily


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
        self._source_index = None
        self._annotation_index = None
        self._system_getter = {}
        self.name = name
        self.taxid = taxid

    def add_file(self, pangenome_file):
        """Links an HDF5 file to the pan. If needed elements will be loaded from this file,
        and anything that is computed will be saved to this file when
        :func:`ppanggolin.formats.writeBinaries.writePangenome` is called.
        :param pangenome_file: A string representing the filepath to the hdf5 pan file
        to be either used or created
        :type pangenome_file: str
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
    def annotation_source(self) -> Set[str]:
        """returns all the annotation source in the pangenomes

        :return: set of annotation source
        """
        source_set = set()
        for gf in self.gene_families:
            for source_annotation in gf.sources:
                source_set.add(source_annotation)
        return source_set

    @property
    def annotations(self) -> Generator[Annotation, None, None]:
        """Create a generator with all annotations in the pangenome

        :return: set of annotation source
        """
        for gf in self.gene_families:
            yield gf.annotations

    def get_gf_by_annnotation(self, annotation: str = None, accession: str = None) -> Generator[GeneFamily, None, None]:
        """ Get gene famlies with a specific annotation or source in pangenome

        :param annotation: Name of the annotation
        :param accession: Accesion identifier of the annotation

        :return: Gene families with the annotation or source
        """
        assert annotation is not None and accession is not None

        for fam in self.gene_families:
            if len(list(fam.get_annotations(name=annotation, accession=accession))) > 0:
                yield fam

    def get_gf_by_sources(self, source: List[str]):
        """ Get gene famlies with a specific source in pangenome

        :param source: Name of the source

        :return: Gene families with the source
        """
        for fam in self.gene_families:
            if fam.get_source(source) is not None:
                yield fam

    @property
    def systems(self) -> Generator[System, None, None]:
        """Get all systems in pangenome
        """
        for system in self._system_getter.values():
            yield system

    def add_system(self, system: System):
        """Add a detected system in the pangenome

        :param system: Detected system that will be added
        """
        same_sys = False
        for system_in in self.systems:
            if system_in.name == system.name or system.name in system_in.canonical_models():
                if system_in.gene_families.issubset(system.gene_families):
                    system.ID = system_in.ID
                    system.add_canonical(system_in)
                    self._system_getter[system.ID] = system
                    same_sys = True

        if not same_sys:
            system.ID = len(self._system_getter) + 1
            self._system_getter[system.ID] = system

    def number_of_systems(self) -> int:
        """Get the number of systems in the pangenomes"""
        return len(self._system_getter)


class Pangenomes:
    """
    This class represente a group of pangenome object.
    """

    def __init__(self):
        """Constructor method
        """
        self.pangenomes_set = {}

    def add_pangenome(self, pangenome: Pangenome):
        """ Add a pangenome object

        :param pangenome: Pangenome object
        :return:
        """
        self.pangenomes_set[pangenome.name] = pangenome
