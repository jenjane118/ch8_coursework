#!/usr/bin python3


""" Gene translation Program """

"""
Program:        getTranslation
File:           getTranslation.py

Version:        1.0
Date:           20.04.18
Function:       Translate coding sequence of gene into peptide sequence.

Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
This program returns peptide translation for specified gene. 


Usage:
======

Revision History:
=================
V1.0           20.04.18             Original            By:JJS
V1.1           23.04.18             Debuggin               JJS

"""
#*****************************************************************************
# Import libraries

import seq_module


#****************************************************************************

# dummy data to try
exon_list = [('AB12345', 36, 807), ('AB12345', 2010, 2222)]
acc = 'AB12345'

with open('seq_file1.txt', 'r') as f:
    file = f.read().splitlines()

#format sequence string
for s in file:
    file = s.replace(' ', '')
    file = s.upper()

## use codingSeq function to get coding sequence
code_seq = seq_module.codingSeq(acc, file, exon_list)

## use tranlate function to get codon list and amino acid translation
peptide = seq_module.translate(acc, code_seq)

##print with codons lined up with amino acids for error checking
for x in peptide[0]:
   print(x, end=' ')
print('')
for x in peptide[1]:
   print( x, end='   ')

# write results to file in xml format

xml_translate = """
    <gene>
        <sequence acc=%(acc)s>
            <codon_seq>%(codon_seq)s</codon_seq>
            <amino_acid_seq>%(amino_acid_seq)s</amino_acid_seq>
        </sequence>
    </gene>
    """
translate_map = {'acc': acc, 'codon_seq': peptide[0], 'amino_acid_seq': peptide[1]}

with open('translate_output.xml', 'w') as xml_file:
    print('<?xml version="1.0" encoding="UTF-8"?>', file=xml_file)
    print(xml_translate % translate_map, file=xml_file)
xml_file.close()
