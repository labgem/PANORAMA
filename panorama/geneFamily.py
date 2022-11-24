#!/usr/bin/env python3
# coding: utf8

# default libraries

# installed libraries
from ppanggolin.geneFamily import GeneFamily as Fam

# local libraries
from panorama.annotation.annotation import Annotation


class GeneFamily(Fam):
    """
    This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.
    :param family_id: The internal identifier to give to the gene family
    :param name: The name of the gene family (to be printed in output files)
    """

    def __init__(self, family_id: int, name: str):
        super().__init__(family_id, name)
        self.annotation = {}  # Key = source, Value = ordered list of best annotation for one source

    def __repr__(self):
        return f"GF {self.ID}: {self.name}"

    def get_annot(self, source: str):
        if self.annotation.get(source) is not None:
            return self.annotation[source]

    def add_annotation(self, source: str, annotation: Annotation, max_prediction: int = None):
        """ Add annotation to gene family

        :param source: Name of database source
        :param annotation: Identifier of the annotation
        """
        source_annot = self.get_annot(source)
        if source_annot is not None:
            index_annot = 0
            insert = False
            while index_annot < len(source_annot):
                current_annot = source_annot[index_annot]
                if current_annot.score is not None and annotation.score is not None:
                    if current_annot.score < annotation.score:
                        source_annot.insert(index_annot, annotation)
                        insert = True
                    elif current_annot.score == annotation.score:
                        if current_annot.e_val is not None and annotation.e_val is not None:
                            if current_annot.e_val < annotation.e_val:
                                source_annot.insert(index_annot, annotation)
                                insert = True
                elif current_annot.e_val is not None and annotation.e_val is not None:
                    if current_annot.e_val < annotation.e_val:
                        source_annot.insert(index_annot, annotation)
                        insert = True
                if not insert:
                    index_annot += 1
                else:
                    break
            if not insert and len(source_annot) < max_prediction:
                source_annot.append(annotation)
        else:
            self.annotation[source] = [annotation]
