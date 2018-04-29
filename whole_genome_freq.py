
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

import gene_module
import seq_module
import codon_usage

from xml.dom import minidom

#****************************************************************************

gene = 'AB12345'


total_freq = {}
## once we have data from database:
#for object in gene_module.Gene._registry:
## use coding seq function to determine coding sequence for each gene
coding_dna = seq_module.codingSeq(gene)

##  call function to determine codon frequency for each gene
codon_table = {}
codon_table = codon_usage.codonFreq(gene, coding_dna)
print(codon_table)

##  add each to total codon frequency dictionary
for key in codon_table:
    if key in total_freq:
        total_freq[key] += codon_table[key]
    else:
        total_freq[key] = codon_table[key]

# calculate codon usage ratio for whole genome (returns dictionary, 'whole_genome_ratio')
gene = 'total'

whole_genome_ratio = codon_usage.usageRatio(gene, total_freq)
print(whole_genome_ratio)

## calculate codon usage percent (usage per 100bp) for whole genome (returns dictionary, 'whole_genome_percent')

whole_genome_percent = codon_usage.codonPercent(gene, total_freq)
print(whole_genome_percent)


#write to file

# freq_xml = """
# <genome>
#     <codon_usage>
#         <codon_freq>%(total_freq)s</codon_freq>
#         <codon_ratio>%(whole_genome_ratio)s</codon_ratio>
#         <codon_percent>%(whole_genome_percent)s</codon_percent>
#     </codon_usage>
# </genome>
# """
# whole_genome_map = {'total_freq':total_freq, 'whole_genome_ratio':whole_genome_ratio, 'whole_genome_percent':whole_genome_percent}

f = open('whole_genome_usage.txt', 'w')
#print('<?xml version="1.0" encoding="UTF-8"?>', file=f)
#print(freq_xml % whole_genome_map, file=f)

print('\n', whole_genome_ratio, file=f)
print('\n', whole_genome_percent, file=f)
f.close()