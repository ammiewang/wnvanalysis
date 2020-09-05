#!/usr/bin/env pypy
"""
Uses the x ecotypes demarcated by ES2 to determine what <x or >x ecotypes would look like using the same sequences.
In the case of < x ecotypes, the back() function is used to walk backwards from the ES2 ecotypes and
merge the most newly divergent clades. In the case of >x ecotypes, the fwd() function is used to walk forwards
on the tree and split the clades closest to the root.
"""
from ete3 import Tree
from Bio import SeqIO
import logparser
from ecotypes import Ecotypes

def set_ecs(ecs, t):
    """creates pairs of (ecotype, common ancestor of ecotype)"""
    s = []
    i = 0
    while i < len(ecs):
        s.append((ecs[i], t.get_common_ancestor(ecs[i])),)
        i += 1
    return s

def get_names(t):
    """get the IDs of all sequences in a tree"""
    new_node = []
    if not t.is_leaf():
        for clade in t.traverse():
            if clade.name:
                new_node.append(clade.name)
    else:
        new_node.append(t.name)
    return new_node


def back(fas, nwk, ys, t):
    """outputs less ecotypes than the original, moves backwards on tree"""
    j = len(ys)
    num_ecotypes = int(input("number of desired ecotypes: "))
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
        keymax = max(temp, key= lambda x: x[2]) #find the putative ecotype in temp whith the greatest distance from the root
        #(i.e. the most newly divergent putative ecotype), we will take this as a new ecotype and repeat the above process until num_ecotypes is reached
        to_remove = []
        m = 0
        while m < len(ys): #find and remove any child nodes of the new ecotype from the ecotype list
            if set(ys[m][0]).issubset(keymax[0]):
                to_remove.append(ys[m])
            m += 1
        for r in to_remove:
            ys.remove(r)
        ys.append((list(keymax[0]),keymax[1]),) #removed ecotypes are merged and appended to the ecotype list
        j = len(ys)
    less_ecotypes = Ecotypes(fas, nwk)
    less_ecotypes.ecotypes = [y[0] for y in ys]
    return less_ecotypes

def fwd(fas, nwk, ys, t):
    """outputs more ecotypes than the original, moves forwards on tree"""
    j = len(ys)
    zs = []
    k = 0
    num_ecotypes = int(input("number of desired ecotypes: "))
    while k < len(ys): #current distance to each ecotype found
        zs.append((ys[k][0], ys[k][1], t.get_distance(ys[k][1])),)
        k += 1
    while j < num_ecotypes:
        xs = []
        for z in zs:
            if len(z[1].children) >= 2 and len(z[0]) > 1: #find all nodes/ecotypes that can be separated into child nodes/ecotypes
                xs.append(z)
        keymin = min(xs, key= lambda x: x[2]) #find the node closest to the root and split it
        ch = keymin[1].children
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
    more_ecotypes = Ecotypes(fas, nwk)
    more_ecotypes.ecotypes = [z[0] for z in zs]
    return more_ecotypes

def to_log(ecotypes):
    ecotypes.original_txt = input("name for log file: ")
    new_ecotypes = open(ecotypes.original_txt, "w")
    count = 1
    ecotype_list = []

    for ecotype in ecotypes.ecotypes:
        new_ecotypes.write("Ecotype " + str(count) + ": " + "[" + ", ".join(x for x in ecotype) + "]" + "\n")
        count += 1
