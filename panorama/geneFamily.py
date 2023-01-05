#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Generator, List, Tuple, Union

# installed libraries
from ppanggolin.geneFamily import GeneFamily as Fam

# local libraries
from annotation import Annotation


class GeneFamily(Fam):
    """
    This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.
    :param family_id: The internal identifier to give to the gene family
    :param name: The name of the gene family (to be printed in output files)
    """

    def __init__(self, family_id: int, name: str):
        super().__init__(family_id, name)
        self._annotationGetter = {}  # Key = source, Value = ordered list of the best annotation for one source
        self.hmm = None
        self.profile = None
        self.optimized_profile = None

    def __repr__(self):
        return f"GF {self.ID}: {self.name}"

    @property
    def annotations(self) -> Generator[Annotation, None, None]:
        """Generate annotations in gene families

        :return: Gene family annotation"""
        for annot_list in self._annotationGetter.values():
            for annotation in annot_list:
                yield annotation

    @property
    def sources(self) -> List[str]:
        """ Get all source annotation in gene family

        :return: List of annotation source
        """
        return list(self._annotationGetter.keys())

    def max_annotation_by_source(self) -> Tuple[str, int]:
        """Get the maximum number of annotation for one source

        :return: Name of the source with the maximum annotation and the number of annotation corresponding
        """
        max_annot = 0
        max_source = None
        for source, annotations in self._annotationGetter.items():
            if len(annotations) > max_annot:
                max_annot = len(annotations)
                max_source = source
        return max_source, max_annot

    def get_source(self, name: str) -> Union[List[Annotation], None]:
        """ Get the annotation for a specific source in gene family

        :param name: Name of the source

        :return: All the annotation from the source if exist else None
        """
        return self._annotationGetter[name] if name in self.sources else None

    def get_annotations(self, name: Union[List[str], str], accession: Union[List[str], str]) -> Generator[Annotation, None, None]:
        """Get annotation by name or accession in gene family

        :param name: Names of annotation searched
        :param accession: Accession number of annotation searched

        :return: annotation searched
        """
        assert name is not None and accession is not None
        name = name if isinstance(name, list) else [name]
        accession = accession if isinstance(accession, list) else [accession]

        for annotation in self.annotations:
            if annotation.name in name or annotation.accession in accession:
                yield annotation

    def add_annotation(self, source: str, annotation: Annotation, max_prediction: int = None):
        """ Add annotation to gene family

        :param source: Name of database source
        :param annotation: Identifier of the annotation
        :param max_prediction:
        """
        for annot in self.get_annotations(name=annotation.name, accession=annotation.accession):
            if annot.source == source:
                self.get_source(source).remove(annot)

        source_annot = self.get_source(source)
        if source_annot is not None:
            index_annot = 0
            insert_bool = False
            while index_annot < len(source_annot):
                current_annot = source_annot[index_annot]
                if current_annot.score is not None and annotation.score is not None:
                    if current_annot.score < annotation.score:
                        source_annot.insert(index_annot, annotation)
                        insert_bool = True
                    elif current_annot.score == annotation.score:
                        if current_annot.e_val is not None and annotation.e_val is not None:
                            if current_annot.e_val > annotation.e_val:
                                source_annot.insert(index_annot, annotation)
                                insert_bool = True
                elif current_annot.e_val is not None and annotation.e_val is not None:
                    if current_annot.e_val > annotation.e_val:
                        source_annot.insert(index_annot, annotation)
                        insert_bool = True
                if not insert_bool:
                    index_annot += 1
                else:
                    break
            if not insert_bool:
                if max_prediction is None or len(source_annot) < max_prediction:
                    source_annot.append(annotation)
            else:
                if max_prediction is not None and len(source_annot) > max_prediction:
                    del source_annot[max_prediction:]
        else:
            self._annotationGetter[source] = [annotation]
