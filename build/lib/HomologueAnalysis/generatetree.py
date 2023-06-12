from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

from retrieve_gene_list import find_homologues

#get list of sequences from genes

if __name__ == "__main__":

    name, animal, sequence = find_homologues("A1BG")

    print(name, animal, sequence)