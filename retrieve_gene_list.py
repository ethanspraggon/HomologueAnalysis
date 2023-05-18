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

#from ensembl import ensembl_coevolution
from ensembl import homologues
from ensemble_query import query

f = open('data.json')

data = json.load(f)

df = pd.DataFrame(data['response']['docs'])

df['Unique Species Names'] = np.nan

df['Sequence Count'] = np.nan





def func(gene_name):
    seq_rec, animal_species, query_name = homologues(gene = value, idcut=0.0)
    return animal_species

animal_dict = {}  
for index, value in df['symbol'].items():
  print("index", index, "value:" ,value)
  print(type(value))
  #continue

  query_list = [value]
  seq_rec, animal_species, query_name = homologues(gene = query_list, idcut=0.0)
  print("seq_req", type(seq_rec))


  #loop through the sequence records, output description, scrape the percent?

  print("desc", type(seq_rec[0].description))

  print("animal_species", animal_species) 
  print("query_name", query_name)

  #create a dictionary with a list of animal species for each gene correlated
  animal_dict[value] = animal_species

  

  percent_matrix_dict[value] = 

  


  #df.at[index, 'Unique Species Names'] = animal_species

  #df.at[index, 'Sequence Count'] = len(seq_rec)

  break


animal_genes_df = pd.DataFrame(animal_dict)

print(animal_genes_df)
  

#percent_matrix_df = pd.DataFrame(columns=[df['symbol']], index=)




#animal_dict = {gene:animal_species for idx, gene in enumerate(df["symbol"])}

  




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