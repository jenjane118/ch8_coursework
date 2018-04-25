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

for s in file:
    file = s.replace(' ', '')
    file = s.upper()

## dummy data
exon_list = [('AB12345.1', 36, 80), ('AB12345.1.1', 84, 95)]
acc = 'AB12345.1'

## get annotated sequence with exon boundaries indicated
ann_seq = seq_module.annotateSeq(acc, file, exon_list)

## get coding sequence
code_seq = seq_module.codingSeq(acc, file, exon_list)

## get translation
aa_seq = seq_module.translate(acc, code_seq)

xml_seq_ann = """
<gene>
    <sequence acc=%(acc)s>
        <ann_seq>%(sequence)s</ann_seq>
    </sequence>
</gene>
"""
seq_map = {'acc': acc, 'sequence': ann_seq}

xml_seq_coding = """
<gene>
    <sequence acc=%(acc)s>
        <coding_seq>%(code_seq)s</coding_seq>
    </sequence>
</gene>
"""
seq_map2 = {'acc':acc, 'code_seq': code_seq}

xml_seq_translation = """
<gene>
    <sequence acc=%(acc)s>
        <aa_seq>%(aa_seq)s</aa_seq>
    </sequence>
</gene>
"""
seq_map3 = {'acc':acc, 'aa_seq':aa_seq[1]}


## print out sequence divided into codons (for lining up with translation)
xml_seq_codons = """
<gene>
    <sequence acc=%(acc)s>
        <codon_seq>%(codon_seq)s</codon_seq>
    </sequence>
 </gene>
"""
seq_map4 = {'acc':acc, 'codon_seq':aa_seq[0]}



seq_map = {'acc': acc, 'sequence': ann_seq}
with open('getSequence_out.xml', 'w') as xml_file:
    print('<?xml version="1.0" encoding="UTF-8"?>', file=xml_file)
    print(xml_seq_ann % seq_map, file=xml_file)
    print(xml_seq_coding % seq_map2, file=xml_file)
    print(xml_seq_translation % seq_map3, file=xml_file)
    print(xml_seq_codons % seq_map4, file=xml_file)
xml_file.close()
