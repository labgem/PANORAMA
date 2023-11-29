#!/usr/bin/env python3
# coding: utf8
import logging

# default libraries

# installed libraries
from ppanggolin.geneFamily import GeneFamily as Fam
from pyhmmer.plan7 import HMM


# local libraries


class GeneFamily(Fam):
    """
    This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.

    :param family_id: The internal identifier to give to the gene family
    :param name: The name of the gene family (to be printed in output files)
    """

    def __init__(self, family_id: int = None, name: str = None, family: Fam = None):
        if family is not None:  # Allow to cast
            if not isinstance(family, Fam):
                raise TypeError("PPanGGOLiN Gene Family type is expected")
            if family_id is not None or name is not None:
                logging.getLogger("Panorama").warning("The provided gene family identifier and name will not be used")
            super().__init__(family.ID, family.name)
            self._edges = family.edges
            self._genePerOrg = family._genePerOrg
            self._genes_getter = family._genes_getter
            self.removed = family.removed  # for the repeated family not added in the main graph
            self.sequence = family.sequence
            self.partition = family.partition
            self._spots = family.spots
            self._modules = family.modules
            self.bitarray = family.bitarray
        elif family_id is not None and name is not None:
            super().__init__(family_id, name)
        else:
            raise AssertionError("You must provide a Gene Family or a name and a family ID")
        self._hmm = None
        self.profile = None
        self.optimized_profile = None

    def __repr__(self):
        return f"GF {self.ID}: {self.name}"

    @property
    def HMM(self) -> HMM:
        """Return gf HMM"""
        return self._hmm

    @HMM.setter
    def HMM(self, hmm: HMM):
        if not isinstance(hmm, HMM):
            raise TypeError(f"Expected type is {HMM.__class__.name}, found type was {type(hmm)}")
        self._hmm = hmm


    # def add_annotation(self, source: str, annotation: Annotation, max_prediction: int = None):
    #     """ Add annotation to gene family
    #
    #     :param source: Name of database source
    #     :param annotation: Identifier of the annotation
    #     :param max_prediction:
    #     """
    #     source_annot = self.get_source(source)
    #     same_name = False
    #     if source_annot is not None:
    #         index_annot = 0
    #         insert_bool = False
    #         while index_annot < len(source_annot):
    #             current_annot = source_annot[index_annot]
    #             if current_annot.name == annotation.name:
    #                 same_name = True
    #             if current_annot.score is not None and annotation.score is not None:
    #                 if current_annot.score < annotation.score:
    #                     if same_name:
    #                         source_annot[index_annot] = annotation
    #                     else:
    #                         source_annot.insert(index_annot, annotation)
    #                         insert_bool = True
    #                 elif current_annot.score == annotation.score:
    #                     if current_annot.e_val is not None and annotation.e_val is not None:
    #                         if current_annot.e_val > annotation.e_val:
    #                             if same_name:
    #                                 source_annot[index_annot] = annotation
    #                             else:
    #                                 source_annot.insert(index_annot, annotation)
    #                                 insert_bool = True
    #             elif current_annot.e_val is not None and annotation.e_val is not None:
    #                 if current_annot.e_val > annotation.e_val:
    #                     if same_name:
    #                         source_annot[index_annot] = annotation
    #                     else:
    #                         source_annot.insert(index_annot, annotation)
    #                         insert_bool = True
    #             if not insert_bool and not same_name:
    #                 index_annot += 1
    #             else:
    #                 break
    #         if not insert_bool:
    #             if not same_name:
    #                 if max_prediction is None or len(source_annot) < max_prediction:
    #                     source_annot.append(annotation)
    #         else:
    #             if max_prediction is not None and len(source_annot) > max_prediction:
    #                 del source_annot[max_prediction:]
    #     else:
    #         self._annotationGetter[source] = [annotation]
