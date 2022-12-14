#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Dict, Generator, List, Union

# install libraries
from ppanggolin.pangenome import Pangenome as Pan

# local libraries
from panorama.annotation import Annotation, System
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

    def get_gene_family(self, name: str) -> GeneFamily:
        return super().get_gene_family(name)
    # def get_gene_family(self, name: str) -> Union[GeneFamily, None]:
    #     try:
    #         fam = super().get_gene_family(name)
    #     except KeyError:
    #         return None
    #     else:
    #         return fam

    @property
    def annotation_source(self) -> set:
        source_set = set()
        for gf in self.gene_families:
            for source_annotation in gf.sources:
                source_set.add(source_annotation)
        return source_set

    @property
    def annotations(self) -> Generator[Annotation, None, None]:
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

    def get_gf_by_sources(self, annotation: str = None, source: List[str] = None):
        """ Get gene famlies with a specific annotation or source in pangenome

        :param annotation: Name of the annotation
        :param source: Name of the source

        :return: Gene families with the annotation or source
        """
        assert source is not None

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
        same_sys = False

        for system_in in self.systems:
            if system_in.name == system.name or system.name in system_in.canonical:
                if system_in.gene_families.issubset(system.gene_families):
                    system.ID = system_in.ID
                    self._system_getter[system_in.ID] = system
                    same_sys = True

        if not same_sys:
            system.ID = len(self._system_getter) + 1
            self._system_getter[system.ID] = system





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
