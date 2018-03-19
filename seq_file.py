
#!/usr/bin python3

"""Return Sequence"""

"""
Program:    seq_file
File:       seq_file.py

Version:    1.0
Date:       19.03.18
Function:   Retrieve sequence for a given accession number

Author:     Jennifer J. Stiens

Course:     MSc Bioinformatics, Birkbeck University of London
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

with open('seq_file.txt', 'r', encoding='utf-8') as f:
    file = f.read().splitlines()

seq_list = []
seq = {}
for line in file:
    acc = line[:7]
    dna = line[9:]
    seq_list.append(dna)
        
with open ('seq_file_out.txt', 'w') as text_file:
    print('<?xml version="1.0" encoding="UTF-8"?>'), file=text_file)
 
    for x in seq_list:
        seq = {'sequence': x}
        print(xml_seq%seq, file=text_file)

    



