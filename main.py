from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
import random

from ecotypes import Ecotypes
from draw_lines import set_ecs, get_names, back, fwd, to_log

if __name__ == '__main__':

#1. draw_lines -----------------------------------------------------------------------------

    print("DRAWING LINES")

    fasta_file = input("fasta file: ")
    nwk_file = input("newick file: ")
    txt_file = input('text file: ')
    orig_ecs = Ecotypes(fasta_file, nwk_file, txt_file)
    orig_ecs.nwk.convert_to_ultrametric()
    orig_ecs.nwk.write(format=1, outfile="ultrametric.nwk")

    inp = input("Would you like less or more ecotypes? Enter 0 for less, 1 for more: ")
    s = set_ecs(orig_ecs.ecotypes, orig_ecs.nwk)
    first_time = True

    while (inp != '0' and inp != '1') or first_time:
        if inp == '0':
            new_ecotypes = back(fasta_file, nwk_file, s, orig_ecs.nwk)
            break
        elif inp == '1':
            new_ecotypes = fwd(fasta_file, nwk_file, s, orig_ecs.nwk)
            break
        else:
            inp = input("Would you like less or more ecotypes? Enter 0 for less, 1 for more: ")
        first_time = False

    to_log(new_ecotypes)

#2. ecotype_files -------------------------------------------------------------------------

    print('\n' + "MAKING ECOTYPE FASTA FILES")

    new_ecotypes.line()

#3. pd_overall -----------------------------------------------------------------------------

    print('\n' + 'CALCULATING EFFECTIVE POPULATION SIZE')

    print("Average Pairwise Divergence in Each Ecotype: ")
    average_pd_across_ecotypes = new_ecotypes.big_file()
    print("Average Pairwise Divergence Across All Ecotypes: ")
    print(average_pd_across_ecotypes)
    print("Effective Population Size: ")
    print(average_pd_across_ecotypes/(6*10**(-5)))
