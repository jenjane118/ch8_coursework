#!/usr/bin python3


""" Gene translation Program """

"""
Program:        getSequence
File:           getSequence.py

Version:        1.0
Date:           23.04.18
Function:       Retrieve genomic and coding sequences for specified gene.

Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
This program returns genomic sequence with coding regions (exon boundaries) indicated.


Usage:
======


Revision History:
=================
V1.0           23.04.18             Original            By:JJS

"""
#*****************************************************************************
# Import libraries

import sys
import seq_module


#****************************************************************************

## dummy data (sequence in file, acc and exon_list below)
with open ('seq_file1.txt', 'r') as f:
    file = f.read().splitlines()
f.close()

## dummy data
exon_list = [('AB12345.1', 36, 807), ('AB12345.1.1', 854, 950)]
acc = 'AB12345.1'

## get annotated sequence with exon boundaries indicated
ann_seq = seq_module.annotateSeq(acc, file, exon_list)

print(ann_seq)

xml_seq_ann = """
<gene>
    <sequence acc=%(acc)s>
        <ann_seq>%(sequence)s</ann_seq>
    </sequence>
</gene>
"""

seq_map = {'acc': acc, 'sequence': ann_seq}
with open('getSequence_out.xml', 'w') as xml_file:
    print('<?xml version="1.0" encoding="UTF-8"?>', file=xml_file)
    print(xml_seq_ann % seq_map, file=xml_file)
xml_file.close()
