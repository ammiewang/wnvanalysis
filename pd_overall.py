"""
Finds the overall pairwise sequence divergence across all ecotypes.
Takes the weighted average of the average pairwise sequence divergence of each ecotype.
Note: this program runs significantly faster using the pypy compiler
and runs quite slowly using the default python compiler due to the large
number of comparisons made in the pd.py file.
"""
from pd import *
import os
def big_file():
    beginning = input("Input the beginning of the file name (e.g. example for example1.fa, example2.fa, etc.): ")
    extension = input("Input the extension (e.g. .fa, .fas, .fasta): ")
    num_files = len([name for name in os.listdir() if name.startswith(beginning) and name.endswith(extension)])
    j = 1
    filenames = []
    while j <= num_files:
        filenames.append(beginning + str(j) + extension) #gathers all the files for processing
        j += 1
    over_one = 0
    sum = 0
    for fname in filenames:
        x, e_len = pd(fname) #submits the ecotypes group by group to the pairwise sequence divergence calculator, receives the PD and the size of the ecotype
        if x is not None: #makes sure to exclude any singleton ecotypes
            sum += x*e_len #adds the weighted PD to the total PD across all ecotypes
            over_one += e_len #if the ecotype is a non-singleton its size is added to the total size count
    return sum/over_one #this returns the weighted average pairwise sequence divergence across all ecotypes

print("Average Pairwise Divergence in Each Ecotype: ")
average_pd_across_ecotypes = big_file()
print("Average Pairwise Divergence Across All Ecotypes: ")
print(average_pd_across_ecotypes)
print("Effective Population Size: ")
print(average_pd_across_ecotypes/(6*10**(-5)))
