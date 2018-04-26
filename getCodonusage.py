
# !/usr/bin python3

""" Codon Usage program """

"""
Program:        getCodonusage
File:           getCodonusage.py

Version:        1.0
Date:           26.04.18
Function:       Returns codon frequency, codon usage ratio and percentage for a particular gene.

Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
This program will calculate codon usage frequency, percentage, and ratio of codon usage preference for a particular gene
given the accession number. It returns 3 dictionaries.


Usage:
======


Revision History:
=================

V1.0           26.04.18         Original                                By: JJS

"""
# *****************************************************************************
# Import libraries

import sys
import seq_module
import codon_usage

# ****************************************************************************


gene = 'AB12345.5'

## obtain coding sequence for gene indicated

code_seq = seq_module.codingSeq(gene)

## calculate raw frequencies of codon usage
codon_freq = codon_usage.codonFreq(gene, code_seq)

## find each amino acid codon usage ratio
ratio       = codon_usage.usageRatio(gene, codon_freq)

## find percent usage (per 100 bp)
percent     = codon_usage.codonPercent(gene, codon_freq)

for k,v in codon_freq.items():
    print(k,v)

for k,v in ratio.items():
    print(k,v)

for k,v in percent.items():
    print(k,v)