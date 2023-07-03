from __future__ import print_function

import argparse
import os
import re
import sys
# third party python modules
###################################################################################################
import pymysql
from Bio import SeqIO
from multiprocessing import Pool, Process
import json
import pickle as pkl
import functools

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#from ensembl import ensembl_coevolution
from HomologueAnalysis.ensembl import homologues
from HomologueAnalysis.ensemble_query import query


from config import cache

#print("cache", cache)
animal_perc_id_dict = {}


def find_homologues(gene_name, id_cut):
    """retrieves homologues from a gene name
    
    Calls the homologues function 
    
    """
    seq_rec, animal_species, query_name = homologues(gene_name, idcut=id_cut)
    #animal_perc_id_dict[gene_name] = animal_species

    return gene_name, animal_species, seq_rec

def generate_percent_matrix(test):

    #exception,
    pool = Pool()
    results = pool.map(functools.partial(find_homologues, id_cut = 0.0), test)

    print("results", results)
    #print("matrix", animal_perc_id_dict)
    return results


if __name__ == "__main__":
    f = open('data.json')

    data = json.load(f)

    df = pd.DataFrame(data['response']['docs'])

    df['Unique Species Names'] = np.nan

    df['Sequence Count'] = np.nan

    


    indexes, names = zip(*list(df['symbol'].items())[:20])

    #get all names that aren't in the cached dictionary

    names = list(names)
    s = set(list(animal_perc_id_dict.keys()))
    new_names = [name for name in names if name not in s]

    results = generate_percent_matrix(names)

    for gene_name, animal_percs, _ in results:
        animal_perc_id_dict[gene_name] = animal_percs

    

    print("matrix", animal_perc_id_dict)

    animal_genes_df = pd.DataFrame(animal_perc_id_dict)

    print(animal_genes_df)

    print("idxmin\n", animal_genes_df.idxmin())

    #generate heatmap, display option
    #phylogenic tree on a gene
    #cut it to vertebrates? print heatmap, set a function to look at a specific gene, highlight it 
    #animal: terrapin_karolina
    #simliarity between dog and human for example, show information how networks compare

    #plot_heatmap(animal_genes_df)

    #red to blue
    animal_genes_df.to_pickle("./animal_gene_matrix.pkl")


    #sns.heatmap(animal_genes_df, cmap=sns.color_palette("viridis", as_cmap=True))

    #plt.show()


    #heatmap, select
