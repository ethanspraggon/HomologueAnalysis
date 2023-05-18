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
def coevolution(seq_records, query, verbose=False):
    recs1, specs1 = seq_records[0]
    recs_combined = []
    description = "-".join(query)
    coevolution_file = f"{description}-ensembl.fas"
    nomatch = []
    for spec1 in specs1.keys():
        #print(spec1, recs1[specs1[spec1]].id)
        rec_list = [recs1[specs1[spec1]], ]
        for idx in range(1, len(seq_records)):
            recs2, specs2 = seq_records[idx]
            if spec1 in specs2:
                #print(recs1[specs1[spec1]].id, recs2[specs2[spec1]].id)
                rec_list.append(recs2[specs2[spec1]])
        # combine if the rec_list is the same size as number of queries
        #####################################
        if len(rec_list) == len(seq_records):
            recid = spec1 + "-" + "-".join([rec.id.split("-")[-1] for rec in rec_list])
            seq_rec = SeqRecord(Seq("".join([str(rec.seq) for rec in rec_list]), ), id=recid, description=description)
            recs_combined.append(seq_rec)
            if verbose:
                print(f"ensembl:coevolution {len(recs_combined)} {description} {spec1} matches")
        else:
            if verbose:
                print(f"ensembl:coevolution {description} {spec1} does not match")
            nomatch.append(spec1)
    # write it out
    #########################################
    print(f"ensembl:coevolution {description} {len(recs_combined)} matched {len(nomatch)} did not match out of {len(recs1)}")
    print(f"ensembl:coevolution writing {len(recs_combined)} sequences for {description} to {coevolution_file}")
    with open(coevolution_file, 'w') as handle:
        SeqIO.write(recs_combined, handle, format='fasta')
    # do clustalw perhaps
    ######################################################################
    align_file = f"{Path(coevolution_file).with_suffix('')}.aln"
    clustalw_exe = clustalwexe()
    if not clustalw_exe is None and not Path(align_file).exists():
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile=coevolution_file)
        print(f"ensembl:coevolution running {clustalw_cline}")
        try:
            clustalw_cline()
        except:
            subprocess.run([clustalw_exe, coevolution_file])
        tree_file = f"{Path(coevolution_file).with_suffix('')}.dnd"
        tree = Phylo.read(tree_file, "newick")
        Phylo.draw_ascii(tree)
    # now do matrix
    ######################################################################
    covariance_matrix(align_file, aformat='clustal')

def ensembl_sequence(id="ENSG00000157764"):
    """ensembl_sequence get the sequence coresponding to id

    Example from https://rest.ensembl.org/documentation/info/sequence_id
    """
    ext = f"/sequence/id/{id}?"
    resp = requests.get(ensembl_server+ext, headers={ "Content-Type" : "text/plain"})
 
    if not resp.ok:
        resp.raise_for_status()
        sys.exit()

    #print(resp.text)
    return resp.text

def genelookup(gene, species="homo_sapiens", expand=0, verbose=False):
    ext = f"/lookup/symbol/{species}/{gene}?expand={expand}"
 
    r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "application/json"})
 
    if not r.ok:
        r.raise_for_status()
        sys.exit()
 
    decoded = r.json()
    if verbose:
        print(repr(decoded))
    return decoded


def genetree(gene, species='homo_sapiens'):
    #server = "https://rest.ensembl.org"
    #ext = "/cafe/genetree/id/ENSGT00390000003602?"
    ext = f"/cafe/genetree/member/symbol/{species}/{gene}?"
     
    #r = requests.get(server+ext, headers={ "Content-Type" : "text/x-nh"})
    r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "application/json"})
     
    if not r.ok:
      r.raise_for_status()
      sys.exit()

    #genetree_dic = json.loads(r.text)
    genetree_dic = r.json()
    print(repr(genetree_dic))
    for key, value in genetree_dic.items():
        if isinstance(value, (str, float, int)):
            print(key, value)
        elif isinstance(value, dict):
            for key2, value2 in value.items():
                print(key, key2, type(value2))

def homologues(gene: list, aligned=0, compara='vertebrates', idcut=40.0, qformat='full', species='homo_sapiens', qtype='orthologues', sequence='protein'):
    # internal functions
    #############################################
    def get_seqrec(dic, idcut=40.0, source=False):
        #print(list(dic.keys()))
        seq_rec = None
        spec = dic['species']
        seq = dic['seq'] if 'seq' in dic else dic['align_seq']
        ens_id = dic['protein_id']
        perc_id = dic['perc_id'] if 'perc_id' in dic else 100.0
        if source:
            perc_id = 100.0
        if seq.find("X") >= 0 or perc_id < idcut:
            return seq_rec
        else:
            #print(idx2, spec, seq, ens_id, perc_id, list(dic.keys()))
            name = "-".join([spec, ens_id])
            seq_rec = SeqRecord(Seq(seq, ), description=f"{species}={perc_id:.1f}%", name=spec, id=name)
        return seq_rec
    # do the query gene is tuple
    ###############################################
    gene_name = None
    homologue_dic = None
    for gen in gene:
        ext = f"/homology/symbol/{species}/{gen}?;format={qformat};type={qtype};sequence={sequence};aligned={aligned};compara={compara}"
        r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "application/json"})
        if not r.ok:
            print(f"ensembl:homologues query for {gen} not good perhaps an alternative")
        else:
            homologue_dic = r.json()
            gene_name = gen
            break
    if homologue_dic is None:
        print(f"ensembl:homologues query for {gene} not good returning None")
        return None, None, None
    # extract seq_recs
    #############################################
    got_source = False
    seq_recs = []
    ref_id = "None"
    animals = {}
    for idx, value in enumerate(homologue_dic["data"]):
        print(f"_______________ {gene_name} {idx} _________________")
        if isinstance(value, (str, float, int)):
            print(f"ensembl:homologues  {idx} {value}")
        elif isinstance(value, dict):
            for value2 in value.values():
                if isinstance(value2, list):
                    for idx2, value3 in enumerate(value2):
                        if not got_source:
                            seq_rec = get_seqrec(value3['source'], idcut=idcut, source=True)
                            ref_id = value3['source']['protein_id']
                            if not seq_rec is None:
                                seq_recs.append(seq_rec)
                                animals[species] = len(animals)
                                got_source = True
                        if 'target' in value3:
                            target = value3['target']
                            seq_rec = get_seqrec(target, idcut=idcut)
                            if not seq_rec is None and not seq_rec.name in animals:
                                seq_recs.append(seq_rec)
                                animals[seq_rec.name] = len(animals)
    # finish up
    #############################################
    output_file = f"{gene_name}_{species}_{ref_id}.fas"
    print(f"ensembl:homologues writing {len(seq_recs)} sequences to file {output_file}")
    #with open(output_file, 'w') as handle:
    #    SeqIO.write(seq_recs, handle, format='fasta')
    return seq_recs, animals, gene_name


def ensembl_coevolution():
    # parse the arguments
    ###############################################################################################
    parser = argparse.ArgumentParser(description="Uniprot coevolution query")
    parser.add_argument('-a', '--aligned', action='store_true', help='download sequences as aligned files')
    parser.add_argument('-b', '--breadth', action='store_true', help='if true do all rather than in pairs')
    parser.add_argument('--coquery', type=str, nargs='+', default=[], help='the coevolution queries')
    parser.add_argument('--fasta', default=None, type=str, help="the fasta file")
    parser.add_argument('--file', default=None, type=str, help="the tab delimited file to read")
    parser.add_argument('--gene', action='store_true', help='this is a gene name')
    parser.add_argument('-i', '--identity', default=40.0, type=float, help="the percent identity cutoff")
    parser.add_argument('--minlen', default=0, type=int, help='minimum length of sequence')
    parser.add_argument('--maxlen', default=1000000, type=int, help='maximim length of sequence')
    parser.add_argument('--output', default="test", type=str, help='the output fasta file')
    parser.add_argument('--query', type=str, default=[], nargs='+', help='query strings')
    parser.add_argument('--sequence', default='protein', type=str, help="the sequence type to output")
    parser.add_argument('--species', default='homo_sapiens', type=str, help="the species")
    parser.add_argument(
                        '--type', 
                        default='orthologues', 
                        type=str, 
                        choices=["orthologues", "paralogues"], 
                        help="the search type"
                        )
    parser.add_argument('--verbose', action='store_true', help='get some more verbose output')
    args = parser.parse_args()
    # assemble query all queries are tuples to take care of 
    ######################################################################
    queries = list([tuple((query,)) for query in args.query])
    # TODO different file formats
    #############################################
    if not args.file is None:
        with open(args.file, 'r') as handle:
            for line in handle.readlines():
                if len(line.strip()) > 0:
                    tmp = line.split("\t")
                    query = tuple((tm.split("_")[0].strip() for tm in tmp))
                    # check they are not all the same
                    if len(set(query)) == 1:
                        query = tuple((query[0], ))
                    queries.append(query)
    # assemble
    ######################################################################
    seq_records = []
    complete_queries = []
    for query in queries:
        aligned = 1 if args.aligned else 0
        seq_recs, species, query_name = homologues(query, aligned=aligned, idcut=args.identity, species=args.species, sequence=args.sequence, qtype=args.type)
        if not seq_recs is None:
            print(f"ensembl:ensembl_coevolution number of seq_records {len(seq_recs)} and unique species {len(species)}")
            seq_records.append((seq_recs, species))
            complete_queries.append(query_name)
    return seq_recs, species
    # combine by breadth take first and combine with all others
    ######################################################################
    """if args.breadth:
        flat_queries = []
        for query in queries:
            flat_queries += list(query)
        print(f"____________ {'-'.join(flat_queries)} _________________________")
        coevolution(seq_records, flat_queries, verbose=args.verbose)
    # do in pairs
    ######################################################################
    else:
        seq_record_pairs = [seq_records[0], None]
        query = [complete_queries[0], None]
        for idx in range(1, len(seq_records)):
            seq_record_pairs[1] = seq_records[idx]
            query[1] = complete_queries[idx]
            print(f"____________ {'-'.join(query)} _________________________")
            coevolution(seq_record_pairs, query, verbose=args.verbose)"""


if __name__ == "__main__":
    res = ensembl_sequence()
    print(dir(res))
    print(res.json())
    #print(res.text)
    sys.exit()
    ensembl_coevolution()

