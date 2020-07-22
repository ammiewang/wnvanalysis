from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
import random

import logparser
from draw_lines import set_s, get_names, back, fwd, to_log
from ecotype_files import line
from pd_overall import big_file
from pd import *
import os

if __name__ == '__main__':

#1. draw_lines -----------------------------------------------------------------------------

    print("DRAWING LINES")

    fasta_file = input("fasta file: ")
    sequences = SeqIO.parse(fasta_file, "fasta")
    sequences_list = list(sequences)
    nwk_file = input("newick file: ")
    t = Tree(nwk_file)
    read_length = len(sequences_list[0])
    txt_file = input("log/text file: ")
    ecs = logparser.parser(txt_file)
    t.convert_to_ultrametric()
    t.write(format=1, outfile="ultrametric.nwk")
    root = t.get_tree_root()

    inp = input("Would you like less or more ecotypes? Enter 0 for less, 1 for more: ")
    s = set_s(ecs,t)
    first_time = True
    
    while (inp != '0' and inp != '1') or first_time:
        if inp == '0':
            new_ecotypes = back(s,s,t)
            break
        elif inp == '1':
            new_ecotypes = fwd(s,t)
            break
        else:
            inp = input("Would you like less or more ecotypes? Enter 0 for less, 1 for more: ")
        first_time = False

    to_log(new_ecotypes)

#2. ecotype_files -------------------------------------------------------------------------

    print('\n' + "MAKING ECOTYPE FASTA FILES")

    sequences = SeqIO.parse(fasta_file, "fasta")
    dict = SeqIO.to_dict(sequences)

    line(dict)

#3. pd_overall -----------------------------------------------------------------------------

    print('\n' + 'CALCULATING EFFECTIVE POPULATION SIZE')

    print("Average Pairwise Divergence in Each Ecotype: ")
    average_pd_across_ecotypes = big_file()
    print("Average Pairwise Divergence Across All Ecotypes: ")
    print(average_pd_across_ecotypes)
    print("Effective Population Size: ")
    print(average_pd_across_ecotypes/(6*10**(-5)))
