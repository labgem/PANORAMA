#!/usr/bin/env python3
# coding:utf-8

# installed libraries
from ppanggolin.pangenome import Pangenome


def get(pangenome: Pangenome):
    print(len(pangenome.modules))