from __future__ import print_function

import argparse
import os
import re
import sys
# third party python modules
###################################################################################################
import pymysql
from Bio import SeqIO
# constants
###################################################################################################
tables=["gene","transcript","exon_transcript","exon"]
# functions
###################################################################################################
def Sample_Tables(cursor, limit=5):
    cursor.execute("""show tables""", ())

    results=cursor.fetchall()
    for res in results:
        print("==================================================")
        for key, value in res.items():
            print(key, value)
            query="select * from "+value+" limit 5"
            #print(query)
            cursor.execute(query,())
            samp=cursor.fetchall()
            print(samp)
    return results
    

#list of human gene names?

#


def query(db="homo_sapiens_core_83_38"):
    connection=pymysql.connect(host="ensembldb.ensembl.org", user='anonymous', db=db, \
                            cursorclass=pymysql.cursors.DictCursor)

    print(type(connection))

    cursor=connection.cursor()
    results = Sample_Tables(cursor)
    sys.exit()

    """
    Do a basic gene annotation
    """

    cursor.execute("""select * from gene limit 10""", ())
    results=cursor.fetchall()
    for i, res in enumerate(results):
        #print(i,res)   
        cursor.execute("""select * from transcript where gene_id=%s""",(res["gene_id"]))
        trans=cursor.fetchall()
        print(i,len(trans)) 
    """
    Other important tables will be transcript, exon_transcript exan translation
    gene --> transcript     gene_id, transcript_id
    exon_transcript ---> transcript_id,exon_id

    """


if __name__ == "__main__":
    query()

