# python modules
###################################################################################################
import argparse
from io import StringIO
import json
import multiprocessing as mp
from multiprocessing import Pool, Queue
import os
from pathlib import Path
import subprocess
import sys
from time import sleep
# third party python modules
###################################################################################################
#from Bio.Alphabet import IUPAC depracated
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from Bio import Phylo
from Bio import pairwise2
from Bio.Align.Applications import ClustalwCommandline
import requests
import numpy as np
# first party python modules
###################################################################################################
#from Seqor.seqor_utils import clustalwexe
#from Seqor.seqor_utils import seqor_metadata
#from Seqor.coevolution import covariance_matrix
# constants
###################################################################################################
ensembl_server = "https://rest.ensembl.org"
# functions
###################################################################################################

def homologues(gene: str, 
               aligned=0, 
               compara='vertebrates', 
               idcut=40.0, 
               qformat='full', 
               species='homo_sapiens', 
               qtype='orthologues', 
               sequence='protein', 
               verbose=False,
               outputfile=False
               ):

    parser = argparse.ArgumentParser(description="homologues")



    """homolgues get ensembl homologues

    Query the ensembl database for genes given in the gene parameter. 
    Compare these sequences across different animal species and store the information.
    Args:
        gene        tuple of strings with gene names and alternatives - only processes one
        aligned     output in the alignment fasta format default is not
        compara.    the animals to look at default is vertebrates
        idcut       float the percentage sequence identity to allow
        qformat     query format default is 'full' which format ot output
        species     string the species to compare against
        qtype       string the query type default is 'orthologues could be paralogues
        sequence    string the sequence type default is protein
        verbose     bool be verbose
        outputfile  bool if print to file
    Output:
        seqrecs    list of seqrecs
        animals    dictionary seqrec.name:position in list
        gene_names string the gene names
    TODO: Percent coverage of matches should be done
    """
    # internal functions
    #############################################
    def get_seqrec(dic, idcut=40.0, source=False, verbose=False):
        """get_seqrec function to return the seqrec from the json result

        Args:
            dic         dictionary returned from ensemble query
            idcut       float the percentage sequence identity cutoff to allow
            sources     bool is this the refence protein or not
            verbose     bool be verbose
        Output:
            seq_rec     SeqRecord biopython seq object
            perc_idcut  float the sequence identity
        """
        seq_rec = None
        spec = dic['species']
        seq = dic['seq'] if 'seq' in dic else dic['align_seq']
        ens_id = dic['protein_id']
        if verbose:
            print(f"get_seqrec: dictionary {dic}")
        perc_id = dic['perc_id'] if 'perc_id' in dic else 100.0
        if source:
            perc_id = 100.0
        if seq.find("X") < 0 and perc_id >= idcut:
            name = "-".join([spec, ens_id])
            seq_rec = SeqRecord(Seq(seq, ), description=f"{species}={perc_id:.1f}%", name=spec, id=name)
        return seq_rec, perc_id
    # do the query gene is tuple
    ###############################################
    gene_name = None
    homologue_dic = None

    #include if intended gene is not known and guesses are made
    """
    for gen in gene:
        ext = f"/homology/symbol/{species}/{gen}?;format={qformat};type={qtype};sequence={sequence};aligned={aligned};compara={compara}"
        r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "application/json"})
        if not r.ok:
            print(f"ensembl:homologues query for {gen} not good perhaps an alternative")
        else:
            homologue_dic = r.json()
            gene_name = gen
            break
    """

    #store informations from ensemble query into a dictionary
    ext = f"/homology/symbol/{species}/{gene}?;format={qformat};type={qtype};sequence={sequence};aligned={aligned};compara={compara}"
    r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "application/json"})


    homologue_dic = r.json()
    gene_name = gene
    

    if homologue_dic is None:
        print(f"ensembl:homologues query for {gene} not good returning None")
        return None, None, None
    # extract seq_recs
    #############################################

    got_source = False
    seq_recs = []
    ref_id = ""
    animals = {}

    if "data" in homologue_dic:
        data = homologue_dic["data"][0]
    else:
        return seq_recs, animals, gene_name


    #length of data is 1? no need to do a for loop
    
    print(f"_______________ {gene_name} _________________")
    
    if not isinstance(data, dict):
        return seq_recs, animals, gene_name
    

    

    # each entries of the data['homologies'] list is a dictionary
    # get the source of the first entry, which will be 

    #check that source seq-req has not been stored and check the first first for the source dictionary

    homologies = data['homologies']

    if homologies == []:
        return seq_recs, animals, gene_name

    
    if not got_source and 'source' in homologies[0].keys():
        seq_rec, perc_id = get_seqrec(homologies[0]['source'], idcut=idcut, source=True)
        ref_id = homologies[0]['source']['protein_id']
        if not seq_rec is None:
            seq_recs.append(seq_rec)
            animals[species] = perc_id
            got_source = True

        
    for homologue in homologies:
        if 'target' in homologue.keys():
            seq_rec, perc_id = get_seqrec(homologue['target'], idcut=idcut)
            if not seq_rec is None and not seq_rec.name in animals:
                seq_recs.append(seq_rec)
                animals[seq_rec.name] = perc_id

    # finish up
    #############################################
    if outputfile:
        output_file = f"{gene_name}_{species}_{ref_id}.fas"
        print(f"ensembl:homologues writing {len(seq_recs)} sequences to file {output_file}")
        with open(output_file, 'w') as handle:
            SeqIO.write(seq_recs, handle, format='fasta')
    if verbose:
        for key, value in animals.items():
            idx, perc_id = value
            print(f"ensembl:homologues  {idx} {gene_name} {key} {perc_id:.2f}%")
    return seq_recs, animals, gene_name

if __name__ == "__main__":
    homologues("A1BG")
    homologues("A1BG-AS1")

