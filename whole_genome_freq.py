
#!/usr/bin python3

""" Whole Genome Codon Usage program """

"""
Program:        whole_genome_freq
File:           whole_genome_freq.py

Version:        1.0
Date:           18.04.18
Function:       Returns codon frequency and codon usage ratio and percentage for entire chromosome.

Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
This program will calculate codon usage frequency, percentage, and ratio of codon usage preference for entire chromosome.


Usage:
======


Revision History:
=================


V1.0           18.04.18         Original        By: JJS

"""
#*****************************************************************************
# Import libraries
import sys

import seq_module

import Codon_Usage

import pickle, shelve
#****************************************************************************

# test sequence
#seq1 = 'CTAAGAATGGCGGCGCTGTGTCGGACCCGTGCTGTGGCTGCCGAGAGCCATTTTCTGCGAGTGTTTCTCTTCTTCAGGCCCTTTCGGGGTGTAGGCACTGAGAGTGGATCCGAAAGTGGTAGTTCCAATGCCAAGGAGCCGCGCGCAGGCGGTTTCGCGAGCGCGTTGGAGCGGCACTCGGAGCTTCTACAGAAGGTGGAGCCCCTACAGAAGGGTTCTCCAAAAAATGTGGAATCCTTTGCATCTATGCTGAGACATTCTCCTCTTACACAGATGGGACCTGCAAAGGATAAACTGGTCATTGGACGGATCTTTCATATTGTGGAGAATGATCTGTACATAGATTTTGGTGGAAAGTTTCATTGTGTATGTAGAAGACCAGAAGTGGATGGAGAGAAATACCAGAAAGGAACCAGGGTCCGGTTGCGGCTATTAGATCTTGAACTTACGTCTAGGTTCCTGGGAGCAACAACAGATACAACTGTACTAGAGGCTAATGCAGTTCTCTTGGGAATCCAGGAGAGTAAAGACTCAAGATCGAAAGAAGAACATCATGAAAAATTT'
#gene = 'AB12345'
#codon_freq = codonFreq(gene, seq1)
#print(codon_freq)


# for each gene in gene_list:
    # run seq_module.codingSeq to get coding sequence

# run with dummy data as coding sequence
gene_list = []
with open('seq_file.txt', 'r') as f:
    file = f.read().splitlines()
    for x in file:
        gene_list.append(x)

# extract accession number and sequence from list of genes

codon_table = {}
total_freq  = {}
for x in gene_list:
    gene = str(x[:7])
    seq  = str(x[8:])

    #for s in seq:
    #    seq_format = s.replace(' ', '')

    codon_table = Codon_Usage.codonFreq(acc=gene, dna=seq)
    #print(codon_table)

    #assert isinstance(codon_table, object)
    for key in codon_table:
        if key in total_freq:
            total_freq[key] += codon_table[key]
        else:
            total_freq[key] = codon_table[key]
print(total_freq)

# calculate codon usage ratio for whole genome
gene = 'total'

whole_genome_ratio = Codon_Usage.usageRatio(gene, total_freq)
#for k,v in whole_genome_ratio.items():
#    print(k, ':', v)

whole_genome_percent = Codon_Usage.codonPercent(gene, total_freq)
#for k,v in whole_genome_percent.items():
#    print(k, ':', v)

#write to file (in pickled format)

f = open('whole_genome_usage.dat', 'wb')

pickle.dump(whole_genome_ratio, f)
pickle.dump(whole_genome_percent, f)
f.close()
