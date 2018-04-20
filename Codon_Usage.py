#!/usr/bin python3

""" Codon Usage Module """

"""
Program:        Codon_Usage
File:           Codon_Usage.py

Version:        1.0
Date:           22.03.18
Function:       Returns codon frequency and codon usage ratio and percentage for a particular gene.

Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
This program will calculate codon usage frequency, percentage, and ratio of codon usage preference for a particular gene.


Usage:
======


Revision History:
=================

V1.0           22.03.18         Original                             By: JJS
V1.1           18.04.18         Renamed (from 'Codon_Usage_module')      JJS

"""
#*****************************************************************************
# Import libraries

import sys


#****************************************************************************

def codonFreq(acc, dna):
    """Returns frequency of each codon possibility in particular sequence.
    Input       acc                     accession number
                dna                     coding sequence (dna)

    Output      CodonsDict              dictionary of codon frequencies in sequence
    """

    CodonsDict = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
        'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    # make a list of codon strings
    codon = ''
    codon_list = []
    for x in dna:
        codon += x
        if len(codon) == 3:
            codon_list.append(codon)
            codon = ''

    # Update dictionary for codon frequencies
    for x in codon_list:
         if x in CodonsDict:
             CodonsDict[x] += 1
    return CodonsDict


def codonPercent (acc, freq_table):
    """Returns percentage use of a particular codon in sequence (per 100 bp sequence)
    Input           acc                 Accession number
                    freq_table          Codon frequency dictionary

    Output          percentDict         Dictionary of codon and usage per 100bp
    """

    total_codons = 0
    percentDict = {}
    for codon in freq_table:
        total_codons += freq_table[codon]
    for codon in freq_table:
        percent = (freq_table[codon] / total_codons)*100
        percentDict[codon] = round(percent, 1)
    return percentDict



def usageRatio (acc, freq_table):
    """Returns codon usage ratio.
    Input       acc             Accession number
                freq_table      Codon frequency dict returned from codonFreq function

    Output      aaDict          Dictionary of amino acids: codons used and ratio of codon usage
    """

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

    aa_list = ['C', 'D', 'S', 'Q', 'M', 'N', 'P', 'K', 'T', 'F', 'A', 'G', 'I', 'L', 'H', 'R', 'W', 'V', 'E', 'Y', '_']



    aaDict = {}
    for aa in aa_list:
        if aa in SynCodons:
            codon_list = SynCodons[aa]
            sum = 0
            ratio = 0
            codonDict = {}
            for x in codon_list:
                sum += freq_table[x]
            for x in codon_list:
                ratio = (freq_table[x]) / sum
                codonDict[x] = round(ratio, 2)
        aaDict[aa] = codonDict
    return aaDict

if __name__ == "__main__":
    print("Ran module directly (and did not 'import' it).")


# main

#
    with open ('seq_file1.txt', 'r') as f:
        file = f.read().splitlines()
# #replace spaces in dna sequence
    seq = ''
    for s in file:
        seq = s.replace(' ', '')
    gene = 'AB12345'
    codon_freq = codonFreq(gene, seq)
    print(codon_freq)

##print each amino acid and codon usage ratio
    ratio = usageRatio(gene, codon_freq)
    for k,v in ratio.items():
        print(k, ':', v)

# # print each codon and percent usage (per 100 bp)
    percent = codonPercent(gene, codon_freq)
    for a in percent:
        print(a, percent[a])
    f.close()