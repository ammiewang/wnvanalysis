#!/usr/bin/env pypy
"""
Creates FASTA files for each ecotype which are then processed by pd_overall.
"""
from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
import random
import math
import operator
import numpy
import logparser
fasta_file = input("fasta file: ")
sequences = SeqIO.parse(fasta_file, "fasta")
dict = SeqIO.to_dict(sequences)
sequences_list = list(sequences)

def line():
    txt = input("log/text file: ")
    name = input("name your fastas (do not include .fas extenstion): ")
    ecotypes = logparser.parser(txt)

    i = 1
    for ec in ecotypes:
        xs = []
        for e in ec:
            x = dict.get(e)
            xs.append(x)
        SeqIO.write(xs, name + str(i) + ".fas", "fasta")
        i += 1

if __name__ == "__main__":
    line()
