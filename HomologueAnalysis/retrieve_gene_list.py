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



def find_homologues(gene_name):
    seq_rec, animal_species, query_name = homologues(gene_name, idcut=0.0)
    #animal_perc_id_dict[gene_name] = animal_species

    return gene_name, animal_species, seq_rec

def generate_percent_matrix(test):

    #exception,
    pool = Pool()
    results = pool.map(find_homologues, test)

    
    print("results", results)
    #print("matrix", animal_perc_id_dict)
    return results

#def plot_heatmap(matrix):

    #plt.style.use("seaborn")

    

    #plt.title( "Heatmap Using Seaborn Method")



"""animal_perc_id_dict = {}  
for index, value in df['symbol'].items():
    #print("index", index, "value:" ,value)
    #print(type(value))
    #continue


    #query_list = [value]
    seq_rec, animal_species, query_name = homologues(gene = value, idcut=0.0)
    #print("seq_req", type(seq_rec))

    print(len(animal_species.keys()))
    #loop through the sequence records, output description, scrape the percent?

    #print("desc", type(seq_rec[0].description))

    #print("animal_species", animal_species) 
    #print("query_name", query_name)

    #create a dictionary with a list of animal species for each gene correlated
    animal_perc_id_dict[value] = animal_species



    #percent_matrix_dict[value] = 

    break




  #df.at[index, 'Unique Species Names'] = animal_species

  #df.at[index, 'Sequence Count'] = len(seq_rec)

  


animal_genes_df = pd.DataFrame(animal_perc_id_dict)

print(animal_genes_df)
  

#percent_matrix_df = pd.DataFrame(columns=[df['symbol']], index=)




#animal_dict = {gene:animal_species for idx, gene in enumerate(df["symbol"])}

  
print(animal_genes_df.idxmin())

"""

"""pool = Pool(5)

def PrintNames(name):
  print(name)
  return True

def func(idx):
  return idx

test = range(0, df["symbol"].shape[0])
results = pool.map(func, test)
print(results)

    """

if __name__ == "__main__":
    f = open('data.json')

    data = json.load(f)

    df = pd.DataFrame(data['response']['docs'])

    df['Unique Species Names'] = np.nan

    df['Sequence Count'] = np.nan

    


    indexes, names = zip(*list(df['symbol'].items())[:5])

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
    sns.heatmap(animal_genes_df)

    plt.show()
