#!/usr/bin python3


""" Gene translation module """

"""
Program:        Translate_module
File:           Translate_module.py

Version:        1.0
Date:           06.04.18
Function:       Translate coding sequence of gene into peptide sequence.
                Calculate codon usage frequency for particulary gene.

Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
This program returns peptide translation and codon usage for specified gene.


Usage:
======

Revision History:
=================

"""
#*****************************************************************************
# Import libraries

from seq_module import codingSeq


#****************************************************************************


def translate(acc, dna):
    aa_seq = ''
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
    }
    codon_list = []
    codon = ''
    count = 0
    #for i in range(1, len(dna)):
    for x in dna:
        #count += 1
        codon += x
        if len(codon) == 3:
            #codon_list.append(x[i]:x[i+3])
            codon_list.append(codon)
            codon = ''
        else:
            continue
    for codon in codon_list:
        if codon in codon_table:
            aa_seq += codon_table[codon]
    return aa_seq



exon_list = [('AB12345', 36, 807), ('AB12345', 2010, 2222)]
acc = 'AB12345'

with open ('seq_file1.txt', 'r') as f:
    file = f.read().splitlines()

code_seq = codingSeq(acc, file, exon_list)
print(code_seq)

peptide =  translate(acc, code_seq)

print(peptide)

# will need to use codon start information to know where to start translating in codingSeq function



