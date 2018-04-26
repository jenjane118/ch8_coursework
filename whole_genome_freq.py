
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

import seq_module

import codon_usage

import pickle, shelve
#****************************************************************************

# test sequence
#seq1 = 'CTAAGAATGGCGGCGCTGTGTCGGACCCGTGCTGTGGCTGCCGAGAGCCATTTTCTGCGAGTGTTTCTCTTCTTCAGGCCCTTTCGGGGTGTAGGCACTGAGAGTGGATCCGAAAGTGGTAGTTCCAATGCCAAGGAGCCGCGCGCAGGCGGTTTCGCGAGCGCGTTGGAGCGGCACTCGGAGCTTCTACAGAAGGTGGAGCCCCTACAGAAGGGTTCTCCAAAAAATGTGGAATCCTTTGCATCTATGCTGAGACATTCTCCTCTTACACAGATGGGACCTGCAAAGGATAAACTGGTCATTGGACGGATCTTTCATATTGTGGAGAATGATCTGTACATAGATTTTGGTGGAAAGTTTCATTGTGTATGTAGAAGACCAGAAGTGGATGGAGAGAAATACCAGAAAGGAACCAGGGTCCGGTTGCGGCTATTAGATCTTGAACTTACGTCTAGGTTCCTGGGAGCAACAACAGATACAACTGTACTAGAGGCTAATGCAGTTCTCTTGGGAATCCAGGAGAGTAAAGACTCAAGATCGAAAGAAGAACATCATGAAAAATTT'
#gene = 'AB12345'
#codon_freq = codonFreq(gene, seq1)
#print(codon_freq)

gene_list = []
with open('seq_file.txt', 'r') as f:
    file = f.read().splitlines()
    for x in file:
        gene_list.append(x)


# for each gene in genelist:
    # run seq_module.codingSeq to get coding sequence

## extract accession number and sequence from list of genes
for x in gene_list:
    gene = str(x[:7])
    new_seq  = str(x[9:])
    print(new_seq)

## get exon list for each gene in list
    ## call on file with exon boundaries list for each gene

##  dummy list:
    gene_exons = [('AC34567', 25, 75), ('AC34567', 85, 95)]


## use coding seq function to determine coding sequence for each gene
    for s in new_seq:
        new_seq += s.replace(' ', '')
        new_seq += s.upper()
    coding_dna = seq_module.codingSeq(gene, new_seq, gene_exons)
    print(coding_dna)
##  call function to determine codon frequency for each gene

    codon_table = {}
    codon_table = codon_usage.codonFreq(gene, coding_dna)
#print(codon_table)

##  add each to total codon frequency dictionary
    total_freq = {}
    for key in codon_table:
        if key in total_freq:
            total_freq[key] += codon_table[key]
        else:
            total_freq[key] = codon_table[key]
print(total_freq)

# calculate codon usage ratio/percent for whole genome
gene = 'total'

whole_genome_ratio = codon_usage.usageRatio(gene, total_freq)
for k,v in whole_genome_ratio.items():
    print(k, ':', v)

whole_genome_percent = codon_usage.codonPercent(gene, total_freq)
for k,v in whole_genome_percent.items():
    print(k, ':', v)

#write to file

freq_xml = """
<genome>
    <codon_usage>
        <codon_freq>%(total_freq)s</codon_freq>
        <codon_ratio>%(whole_genome_ratio)s</codon_ratio>
        <codon_percent>%(whole_genome_percent)s</codon_percent>
    </codon_usage>
</genome>
"""
whole_genome_map = {'total_freq':total_freq, 'whole_genome_ratio':whole_genome_ratio, 'whole_genome_percent':whole_genome_percent}

f = open('whole_genome_usage.xml', 'w')
print('<?xml version="1.0" encoding="UTF-8"?>', file=f)
print(freq_xml % whole_genome_map, file=f)
#print('\n', total_freq, file=f)
#print('\n', whole_genome_ratio, file=f)
#print('\n', whole_genome_percent, file=f)
f.close()