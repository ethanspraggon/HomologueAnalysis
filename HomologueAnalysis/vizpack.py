import argparse
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

from retrieve_gene_list import find_homologues
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

#get list of sequences from genes

def generate_heatmap(matrix, wanted_gene, wanted_animal):

    """ Create a heatmap from a percentage identity matrix

    A heatmap is created from a matrix that specifies the animal and genes percentage identity relationship to humans.
    The user can highlight a certain gene and animal to be seen in the heatmap

    args:
        matrix: an matrix of percentage idenity information
        wanted_gene: the gene whose column will be highlighted
        wanted_animal: the animal whose row will be highlighted
    """
    ax = sns.heatmap(matrix, cmap=sns.color_palette("viridis", as_cmap=True), linewidths=.5)

    rows, cols = matrix.shape

    genes = list(matrix.columns)
    animals = matrix.axes[0].tolist()
    
    wanted_animal = wanted_animal
    wanted_gene = wanted_gene

    gene_index = genes.index(wanted_gene)
    animal_index = animals.index(wanted_animal)
    x, y, w, h = gene_index, 0, 1, rows
    ax.add_patch(Rectangle((x,y), w, h, fill=False, edgecolor='crimson', lw=4, clip_on = False))
    ax.tick_params(length = 0)

    x,y,w,h = 0, animal_index, cols, 1

    ax.add_patch(Rectangle((x,y), w, h, fill=False, edgecolor='crimson', lw=0.5, clip_on = False))
    ax.tick_params(length = 0)

    return ax


    

if __name__ == "__main__":

    argparse

    animal_gene_matrix = pd.read_pickle("./animal_gene_matrix.pkl")
    heatmap = generate_heatmap(animal_gene_matrix, "AACS", "naja_naja")


    #name, animal, sequences = find_homologues("A1BG")
    #print(name, animal, sequences)

    #sequences_dict = {sequence.id: str(sequence.seq) for sequence in sequences}
    plt.show()



