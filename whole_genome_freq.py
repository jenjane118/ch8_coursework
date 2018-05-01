
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
whole_genome_freq

Revision History:
=================


V1.0           18.04.18         Original            By: JJS
V1.1           1.05.18          changed output          JJS

"""
#*****************************************************************************
# Import libraries

import gene_module
import seq_module
import codon_usage

from xml.dom import minidom

#****************************************************************************
def help():
    """Print a usage message and exit."""
    print("""
    codon_usage.py   V1.1        2018,   J.J. Stiens

    Usage: 
    ============
    codon_usage 

    Description:
    ============
    Calculate codon usage percentage, and ratio of codon usage preference for entire chromosome
    First dictionary is amino acid: list of synonymous codons
    Second dictionary is codons: ratio, percent
    """)
    exit(0)

#****************************************************************************

total_freq = {}
## once we have data from database:
#for object in gene_module.Gene._registry:
## use coding seq function to determine coding sequence for each gene
## dummy data
gene = 'AB371373.1'

coding_dna = seq_module.codingSeq(gene)

SynCodons = {
        'C': ['TGT', 'TGC'],
        'D': ['GAT', 'GAC'],
        'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'Q': ['CAA', 'CAG'],
        'M': ['ATG'],
        'N': ['AAC', 'AAT'],
        'P': ['CCT', 'CCG', 'CCA', 'CCC'],
        'K': ['AAG', 'AAA'],
        'T': ['ACC', 'ACA', 'ACG', 'ACT'],
        'F': ['TTT', 'TTC'],
        'A': ['GCA', 'GCC', 'GCG', 'GCT'],
        'G': ['GGT', 'GGG', 'GGA', 'GGC'],
        'I': ['ATC', 'ATA', 'ATT'],
        'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
        'H': ['CAT', 'CAC'],
        'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'W': ['TGG'],
        'V': ['GTA', 'GTC', 'GTG', 'GTT'],
        'E': ['GAG', 'GAA'],
        'Y': ['TAT', 'TAC'],
        '_': ['TAG', 'TGA', 'TAA']}

##  call function to determine codon frequency for each gene
codon_table = codon_usage.codonFreq(gene, coding_dna)

##  add each to total codon frequency dictionary
for key in codon_table:
    if key in total_freq:
        total_freq[key] += codon_table[key]
    else:
        total_freq[key] = codon_table[key]

# calculate codon usage ratio for whole genome (returns dictionary, 'whole_genome_ratio')
gene = 'total'

whole_genome_ratio = codon_usage.usageRatio(gene, total_freq)

## calculate codon usage percent (usage per 100bp) for whole genome (returns dictionary, 'whole_genome_percent')
whole_genome_percent = codon_usage.codonPercent(gene, total_freq)

## results as one dictionary with both ratio and percent for each codon:
ratio_list = []
ratio_dict = {}
usage_dict = {}
for k, v in whole_genome_ratio.items():
    ratio_list.append(v)

for item in ratio_list:
    for k,v in item.items():
        ratio_dict[k] = v

for codon, ratio in ratio_dict.items():
    for codon, percent in whole_genome_percent.items():
        usage_dict[codon] = ratio_dict[codon], percent
print(SynCodons)
print(usage_dict)                 #dictionary listing codon: ratio, percent

help()


#write to file
f = open('whole_genome_usage.txt', 'w')
print(SynCodons, file=f)
print(usage_dict, file=f)
f.close()