from ete3 import Tree
from Bio import SeqIO
import os

class Ecotypes:
    def __init__(self, fasta_file, nwk_file, original_txt=None):
        self.fasta_file = fasta_file
        self.nwk_file = nwk_file
        self.sequences = SeqIO.parse(fasta_file, "fasta")
        self.sequences_list = list(self.sequences)
        self.read_length = len(self.sequences_list[0])
        self.nwk = Tree(nwk_file)
        self.root = self.nwk.get_tree_root()
        self.extension = None
        self.ind_name = None
        if original_txt is not None:
            self.original_txt = original_txt
            self.ecotypes = self.parser(original_txt)

    def line(self):
        """
        Creates FASTA files for each ecotype which are then processed by pd_overall.
        """
        self.sequences = SeqIO.parse(self.fasta_file, "fasta")
        self.seq_dict = SeqIO.to_dict(self.sequences)
        self.ind_name = input("name your fastas (do not include .fas extenstion): ")
        self.extension = '.fas'
        i = 1
        for ec in self.ecotypes:
            xs = []
            for e in ec:
                x = self.seq_dict.get(e)
                xs.append(x)
            SeqIO.write(xs, self.ind_name + str(i) + ".fas", "fasta")
            i += 1

    def big_file(self):
        """
        Finds the overall pairwise sequence divergence across all ecotypes.
        Takes the weighted average of the average pairwise sequence divergence of each ecotype.
        Note: this function runs significantly faster using the pypy compiler
        and runs quite slowly using the default python compiler due to the large
        number of comparisons made in the pd function.
        """
        if self.ind_name is None:
            self.ind_name = input("Input the beginning of the file name (e.g. example for example1.fa, example2.fa, etc.): ")
        if self.extension is None:
            self.extension = input("Input the extension (e.g. .fa, .fas, .fasta): ")
        num_files = len([name for name in os.listdir() if name.startswith(self.ind_name) and name.endswith(self.extension)])
        j = 1
        filenames = []
        while j <= num_files:
            filenames.append(self.ind_name + str(j) + self.extension) #gathers all the files for processing
            j += 1
        over_one = 0
        sum = 0
        for fname in filenames:
            x, e_len = self.pd(fname) #submits the ecotypes group by group to the pairwise sequence divergence calculator, receives the PD and the size of the ecotype
            if x is not None: #makes sure to exclude any singleton ecotypes
                sum += x*e_len #adds the weighted PD to the total PD across all ecotypes
                over_one += e_len #if the ecotype is a non-singleton its size is added to the total size count
        return sum/over_one

    def pd(self, f): #pd for pairwise divergence
        """
        Calculates the average pairwise sequence divergence in one ecotype.
        """
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
            return None,  len(sequences_list)

    def parser(self, inp):
        """
        parses a text/log file containing ecotypes into list of ecotypes
        """
        file = open(inp)
        ecotypes = []
        for line in file:
            if "Ecotype" in line and "[" in line and "]" in line:
                start = False
                ecotype = []
                sequence = ""
                for char in line:
                    if char == '[':
                        start = True
                    if start == True:
                        if char != '[' and char != ',' and char != ' ' and char != ']':
                            sequence += char
                        elif char == ' ':
                            ecotype.append(sequence)
                            sequence = ""
                        elif char == ']':
                            ecotype.append(sequence)
                            ecotypes.append(ecotype)
                            sequence = ""
                            break
        return ecotypes
