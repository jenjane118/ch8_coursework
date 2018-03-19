
#!/usr/bin python3

"""Return Sequence"""

"""
Program:    seq_file
File:       seq_file.py

Version:    1.0
Date:       19.03.18
Function:   Retrieve sequence for a given accession number

Author:     Jennifer J. Stiens

Course:     MSc Bioninformatics, Birkbeck University of London
            Biocomputing2 Coursework Assignment

_____________________________________________________________________________

Description:
============
This program returns the raw sequence for a given gene.

Eventually, indicate exon boundaries and restriction enzyme cutting sites

Usage:
======

Revision History:
=================



"""

#******************************************************************************
# Import libraries

import os

#*****************************************************************************

xml_seq = """
<gene>
<sequence>
<rawseq>%(sequence)s</rawseq>
</sequence>
</gene>
"""


## Read text file including accession numbers and sequence, putting in dictionary
## these will probably be in form of tuples(or dictionaries?) returning from pymysql queries

with open('seq_file.txt', 'r', encoding='utf-8') as f:
    file = f.read().splitlines()

seq_list = []
seq = {}
for line in file:
    acc = line[:7]
    dna = line[9:]
    info = (acc, dna)
    seq_list.append(info)
#print(seq_list)

for x in seq_list:
    if x not in seq:
        seq[x[0]] = x[1]
    else:
        pass
    
"""
A function used to retrieve sequence from dictionary of accession numbers and sequences.

    Input:              gene            Accession number input from user interface

    Output:             seq_map         Returns a mapping of sequence for use in xml template

"""

def GeneSeq(gene):
    if gene in seq:
        seq_map = {'sequence' : seq[gene]}
    return seq_map


user_gene = 'AC34569'

sequence = GeneSeq(user_gene)


with open ('seq_file_out.txt', 'w') as text_file:
    print('<?xml version="1.0" encoding="UTF-8"?>', file=text_file)
    print(xml_seq%sequence, file=text_file)

    



