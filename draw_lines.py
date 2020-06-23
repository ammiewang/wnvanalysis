#!/usr/bin/env pypy
"""
Uses the x ecotypes demarcated by ES2 to determine what <x or >x ecotypes would look like using the same sequences.
In the case of < x ecotypes, the back() function is used to walk backwards from the ES2 ecotypes and
merge the most newly divergent clades. In the case of >x ecotypes, the fwd() function is used to walk forwards
on the tree and split the clades closest to the root.
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
sequences_list = list(sequences)
nwk_file = input("newick file: ")
t = Tree(nwk_file)
read_length = len(sequences_list[0])
txt_file = input("log/text file: ")
ecs = logparser.parser(txt_file)
t.convert_to_ultrametric()
t.write(format=1, outfile="ultrametric.nwk")
root = t.get_tree_root()

node = []
for seq in sequences_list:
    node.append(seq.id)

i = 0
ecss = []
while i < len(ecs): #make ecotypes into sets, easier for comparison
    ecss.append(set(ecs[i]))
    i += 1

s =[]
i = 0
while i < len(ecs): #pairs of (ecotype, common ancestor of ecotype)
    s.append((ecs[i], t.get_common_ancestor(ecs[i])),)
    i += 1

def get_names(t): #get the IDs of all sequences in a tree
    new_node = []
    if not t.is_leaf():
        for clade in t.traverse():
            if clade.name:
                new_node.append(clade.name)
    else:
        new_node.append(t.name)
    return new_node


def back(xs,ys,t): #outputs less ecotypes than the original, moves backwards on tree
    j = len(ys)
    num_ecotypes = input("number of desired ecotypes: ")
    while j > num_ecotypes:
        temp = []
        k = 0
        while k < len(ys):
            newt = t.get_common_ancestor(ys[k][0])
            p = newt.up #find next node up from the common ancestor of an ecotype, this will be a new putative ecotype
            if p:
                p_names = get_names(p) #get IDs of all sequences in the node
                dist = t.get_distance(p) #get distance from the root
                if (set(p_names),p,dist) not in temp:
                    temp.append((set(p_names),p,dist),) #add the merged ecotype to temporary list
            k += 1
        Keymax = max(temp, key= lambda x: x[2]) #find the putative ecotype in temp whith the greatest distance from the root
        #(i.e. the most newly divergent putative ecotype), we will take this as a new ecotype and repeat the above process until num_ecotypes is reached
        to_remove = []
        m = 0
        while m < len(ys): #find and remove any child nodes of the new ecotype from the ecotype list
            if set(ys[m][0]).issubset(Keymax[0]):
                to_remove.append(ys[m])
            m += 1
        for r in to_remove:
            ys.remove(r)
        ys.append((list(Keymax[0]),Keymax[1]),) #removed ecotypes are merged and appended to the ecotype list
        j = len(ys)
    return ys

def fwd(ys,t):
    j = len(ys)
    zs = []
    k = 0
    num_ecotypes = input("number of desired ecotypes: ")
    while k < len(ys): #current distance to each ecotype found
        zs.append((ys[k][0], ys[k][1], t.get_distance(ys[k][1])),)
        k += 1
    while j < num_ecotypes:
        xs = []
        for z in zs:
            if len(z[1].children) >= 2 and len(z[0]) > 1: #find all nodes/ecotypes that can be separated into child nodes/ecotypes
                xs.append(z)
        Keymin = min(xs, key= lambda x: x[2]) #find the node closest to the root and split it
        ch = Keymin[1].children
        to_remove = []
        for x in xs:
            for c in ch:
                if set(get_names(c)).issubset(set(x[0])): #remove parent node of the new split nodes
                    if x not in to_remove:
                        to_remove.append(x)
        for r in to_remove:
            zs.remove(r)
        for c in ch:
            zs.append((get_names(c),c,t.get_distance(c)),) #append two new ecotypes to ecotypes
        j += 1
    return [(z[0],z[1]) for z in zs]

def to_log(ecotypes):
    out_file = input("name for log file: ")
    new_ecotypes = open(out_file, "w")
    count = 1
    ecotype_list = []

    for ecotype in ecotypes:
        ecotype_list.append(ecotype[0])

    for ecotype in ecotype_list:
        new_ecotypes.write("Ecotype " + str(count) + ": " + "[" + ", ".join(x for x in ecotype) + "]" + "\n")
        count += 1

if __name__ == "__main__":
    less_ecotypes = back(s,s,t)
    #print(less_ecotypes)
    #print(len(less_ecotypes))
    to_log(less_ecotypes)
    more_ecotypes = fwd(s,t)
    #print(more_ecotypes)
    #print(len(more_ecotypes))
    to_log(more_ecotypes)
