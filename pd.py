#!/usr/bin/env pypy
"""
Calculates the average pairwise sequence divergence in one ecotype.
"""
from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
import random
import math
import operator
import numpy
import logparser

def pd(f): #pd for pairwise divergence
    sequences = SeqIO.parse(f, "fasta")
    sequences_list = list(sequences)
    if len(sequences_list) == 1: #if there is only a single sequence we cannot measure pairwise sequence divergence so we immediately return
        return None, 1

    matrix = [] #the pairwise divergence for sequences i and j will be stored at matrix[i][j] while matrix[j][i] will remain None as they are the same
    for i in range(len(sequences_list)):
        m = []
        for j in range(len(sequences_list)):
            m.append(None)
        matrix.append(m)

    nucleotides = {'A', 'G', 'C', 'T', '-'}
    #this part is very slow, especially if you have large ecotypes
    for i in range(len(sequences_list)):
        ui = sequences_list[i].seq.upper()
        for j in range(i+1, len(sequences_list)):
            x = 0
            seq_len = len(sequences_list[0].seq) #every time we pick a new sequence j to pair with sequence i, we revert to the shared length of all the aligned sequences
            uj = sequences_list[j].seq.upper()
            for k in range(len(sequences_list[0].seq)):
                if ui[k] not in nucleotides or uj[k] not in nucleotides: #the efficiency can be somewhat increased without this check although this may deteriorate accuracy
                    seq_len -= 1 #every time a nucleotide is not in the set {'A', 'G', 'C', 'T', '-'},
                    #we subtract one from the length since we will not consider these nucleotides in the pairwise divergence calculation
                elif ui[k] != uj[k]:
                    x += 1 #every time both nucleotides are in the set {'A', 'G', 'C', 'T', '-'} but do not match, we add one to the "divergence" count
            pd = x / seq_len #we take the divergence count over the length to get the pairwise divergence for sequences i and j
            matrix[i][j] = pd

    sum = 0
    num = 0
    for x in range(len(matrix)):
        for y in range(len(matrix[0])):
            if matrix[x][y] is not None: #across all valid matrix entries we take a sum of the pairwise divergence values and count the number of pairs
                sum += matrix[x][y]
                num += 1
    if num is not 0:
        print(sum/num)
        return sum/num, len(sequences_list) #sum/num is the average pairwise divergence in an ecotype while len(sequences_list) is the size of the ecotype
    else:
        #print(0)
        return None,  len(sequences_list)
