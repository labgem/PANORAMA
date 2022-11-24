#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Union, List

# local libraries


class Annotation:
    """
    This represents an annotation for one or more Gene families.

    :param source: source of the annotation
    :param accession: source accesion identifier
    :param name: name of the annotation/function
    :param secondary_names: Other possible name for the annotation
    :param description: description of the annotation
    :param e_val: E-value of the gene family/profile comparison
    :param score: Bit score of the gene family/profile comparison.
    :param bias: The biased composition score correction that was applied to the bit score
    """

    def __init__(self, source: str, accession: str, name: str,
                 secondary_names: Union[str, List[str]] = None, description: str = None,
                 score: float = None, e_val: float = None, bias: float = None):
        assert any(x is not None for x in [score, e_val])
        self.source = source
        self.accession = accession
        self.name = name
        self.secondary_names = secondary_names
        self.description = description
        self.score = score
        self.e_val = e_val
        self.bias = bias
