#!/usr/bin python3

""" Sequence annotations module """

"""
Program:        Seq_module
File:           Seq_module.py

Version:        1.0
Date:           22.03.18
Function:       Includes all sequence-related functions.
Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
This program retrieves the genomic sequence for a particular gene and returns sequence information,
including:
Annotated sequence with exons indicated
Coding sequence
Peptide translation
Restriction enzyme analysis

Usage:
======
seq_module      

Revision History:
=================
V1.0            22.03.18    Original                        By: JJS
V1.1            06.04.18    Revised as module                   JJS
V1.2            21.04.18    Combined translate function         JJS
V1.2            24.04.18    Debugging of exon annotation        JJS
V1.3            26.04.18    Simplified parameters               JJS
V1.4            01.05.18    Debugging/unkn bp, exon boundaries  JJS
"""
#*****************************************************************************
# Import libraries


import re
import sys

#****************************************************************************

def help():
    """Print usage message and exit"""
    print("""
    seq_module.py           V1.0        2018,            JJStiens
    
    Description:
    ============
    This program retrieves the genomic sequence for a particular gene and returns sequence information,
    including:
    Annotated sequence with exons indicated
    Coding sequence
    Peptide translation
    Restriction enzyme analysis

    Usage:
    ======
    seq_module      

    """)

    exit(0)
#****************************************************************************

def getSequence(acc):
    """ Return genomic DNA sequence in numbered form.
    Input               acc                 Accession ID

    Output              num_seq             Dictionary of numbered bp for requested sequence
    """


    ## use pymysql script with accession number to access database for 'sequence'
    # import sequence

    # gene = sequence[0]
    ## check to see accession number requested matches return
    #if gene == acc:
    #    seq = sequence[1]
    #else:
    #   print('Error: Sequence not found')
     ## exit program?

    ## dummy sequence until have data
    with open ('sequence_test_out.txt', 'r') as f:
        file = f.read().splitlines()
    f.close()
    seq = ''
    for x in file:
        seq += x
    #seq = 'atggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacgagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaa'

    ## format sequence in unbroken, uppercase form (do this in individual programs
    seq = seq.replace(' ', '')
    seq = seq.upper()

    num_seq = {}
    count = 0
    for x in seq:
        count += 1
        num_seq[count] = x
    return num_seq

#**********************************************************************************

def annotateSeq(acc):           ## need to enter sequence, exon_list until get usable data from db
    """Return sequence with exon boundaries marked out with symbols (or simple sequence string if no exons indicated).
    Input               acc                 Accession ID

    Output              exon_seq            Annotated sequence with inserted *exon/exon* boundaries (string)
    """

    ## use pymysql script with accession number to access database for 'sequence'
    # import sequence
    # gene = sequence[0]
    ## check to see accession number requested matches return
    ## extract sequence
    # if gene == acc:
    #    seq = sequence[1]
    # else:
    #   print('Error: Sequence not found')
    ## exit program?

    ## use pymysql script with accession number to retrieve 'exon'
    ## check to see accession number requested matches database return
    ## make list of exons
    # for x in exon_file:
    #    if x[0] = acc:
    #        x[1] = exon
    #        exon_list.append(exon)

    ## dummy sequence until have data
    with open ('sequence_test_out.txt', 'r') as f:
        file = f.read().splitlines()
    f.close()
    seq = ''
    for x in file:
        seq += x
    #seq = 'atggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacgagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaa'

    exon_list = [('AB371373.1', 2125, 2215), ('AB371373.1', 3642, 3728), ('AB371373.1', 6222, 6300),
                 ('AB371373.1', 9012, 9086), ('AB371373.1', 10313, 10358), ('AB371373.1', 11120, 11264)]
    #exon_list = [('U16860.1', 1, 219)]

    exon_seq = ''

    seq = seq.replace(' ', '')
    seq = seq.upper()

    base_count   = 0
    if exon_list    != None:
        for base in seq:
            base_count += 1
            exon_seq += base
            first_exon = True
            for x in exon_list:
                if first_exon == True:
                    start   = x[1]
                    end     = x[2]
                    first_exon = False
                ## correction for exon start after first exon
                else:
                    start   = (x[1]) -1
                    end     = x[2]
                ## insert 'exon' markers at start/end of exons
                if base_count == start:
                    exon_seq += '*exon'
                elif base_count == end:
                    exon_seq += 'exon*'
    else:
        exon_seq = seq
    return exon_seq

#**********************************************************************************

def codingSeq(acc):
    """Return coding sequence (stuck together exons). If no exon_list, will return genomic sequence string.
    Input           acc                 Accession ID

    Output          coding_seq          Coding sequence
    """

    ## use pymysql script with accession number to access database for 'sequence' (make a function?)
    # import sequence
    # gene = sequence[0]
    ## check to see accession number requested matches return
    ## extract sequence
    # if gene == acc:
    #    seq = sequence[1]
    # else:
    #   print('Error: Sequence not found')
    ## exit program?

    ## use pymysqul script with accession number to retrieve 'exon'
    ## make this a function?
    ## check to see accession number requested matches database return
    ## make list of exons
    #for x in exon_file:
    #    if x[0] = acc:
    #        x[1] = exon
    #        exon_list.append(exon)

    ## dummy sequence until have data
    #seq = 'atggcggcgctgtgtcggacccgtgctgggggtaggcgcgcgcgcgcgcggtaggcactgagagtggtacatatttttagggatatcgctcgagagagagagacacacacacacaccacaccacaatatattatattatattagggagaggaatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacgagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaa'
    with open ('sequence_test_out.txt', 'r') as f:
        file = f.read().splitlines()
    f.close()
    seq = ''
    for x in file:
        seq += x

    ## dummy data
    exon_list = [('AB371373.1', 2125, 2215), ('AB371373.1', 3642, 3728), ('AB371373.1', 6222, 6300), ('AB371373.1', 9012, 9086), ('AB371373.1', 10313, 10358), ('AB371373.1', 11120, 11264)]
    #exon_list = [('U16860.1', 1, 219)]
    codon_start = 2

    seq = seq.replace(' ', '')
    seq = seq.upper()

    if exon_list != None:
        coding_seq = ''
        first_exon = True
        # have to subtract two from codon_start to align
        codon_start = codon_start - 2
        if codon_start < 0:
            codon_start = 0
        for x in exon_list:
            if first_exon == True:
                start = x[1]
                start     = start + codon_start
                end         = x[2]
                first_exon = False
            ## correction needed to make exon boundaries work after first exon
            else:
                start   = (x[1]) - 1
                end     = x[2]
            coding_seq += seq[start:end]

    else:
        coding_seq = seq
    return coding_seq

#**********************************************************************************

def translate(acc):
    """Return protein translation for requested gene.
    Input                       acc                     gene accession number

    Output [0]                  codon_list              ordered list of codons
    Output [1]                  aa_seq                  string of amino acid sequence
    """

    ## must run codingSeq function to obtain coding sequence for translation
    seq = codingSeq(acc)

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

    ## divides string into list of codons
    for x in seq:
        codon += x
        if len(codon) == 3:
            codon_list.append(codon)
            codon = ''

    ## finds corresponding amino acid for codon in list
    for codon in codon_list:
        if codon in codon_table:
            aa_seq += codon_table[codon]
        ## if codon includes unknown nucleotides, insert x for unkown amino acid
        else:
            aa_seq += 'x'
    return codon_list, aa_seq

#**********************************************************************************

def enz_cut(acc, seq=None, enzyme=None):
    """ Indicate any cleavage sites from restriction enzymes in
        restriction enzyme dictionary. Optionally, search a custom cleavage site (enzyme).
        Return dictionary of enzyme name and cleavage positions.

    Input           acc                 accession number for gene
                    sequence            gene sequence string (optional, set as 'None' if not using)
                    enzyme              custom cleavage site (optional)

    Output          cut_dict            dictionary {enzyme:no. of cleavage sites, positions}

    """

    enz_dict = {
        'EcoRI': 'GAATTC', 'BamHI': 'GGATCC',
        'BsuMI': 'CTCGAG', 'HindIII': 'AAGCTT',
        'EcoRV': 'GATATC', 'Sma1': 'CCCGGG'}

    ## no sequence parameter: check default, genomic sequence of gene for restriction enzyme sites
    ## sequence parameter will be utilised when checking coding sequence in 'getEnzymes' program
    if seq == None:
        ## use pymysql script with accession number to access database for 'sequence'
        # import sequence
        # gene = sequence[0]
        ## check to see accession number requested matches return
        ## extract sequence
        # if gene == acc:
        #    seq = sequence[1]
        # else:
        #   print('Error: Sequence not found')
        ## exit program?

        ## dummy sequence
        #seq = 'atggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacgagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaa'

        with open('AB007516.1_out.txt', 'r') as f:
            file = f.read().splitlines()
        f.close()
        seq = ''
        for x in file:
            seq += x

        #exon_list = [('AB371373.1', 2125, 2215), ('AB371373.1', 3642, 3728), ('AB371373.1', 6222, 6300),
        #             ('AB371373.1', 9012, 9086), ('AB371373.1', 10313, 10358), ('AB371373.1', 11120, 11264)]

    seq = seq.replace(' ', '')
    seq = seq.upper()

    seq_count = 0
    cut_dict = {}

    ## if a custom cleavage site is included, will indicate number and position of cleavage sites in sequence
    if enzyme != None:
        cut_list = []
        count = 0
        cut = enzyme
        ## find custom cleavage sites if present in sequence
        p = re.compile(r'(' + cut + ')')
        it = p.finditer(seq)

        ## count number of cleavage sites
        for match in it:
            count += 1
            ## if cleavage site is present in sequence, add to cut dictionary
            if count != 0:
                cut_list.append((match.start(), match.end()))
                cut_dict[enzyme] = (count, cut_list)

    ## for each enzyme in dictionary, it will search and return only the enzymes which have cleavage sites in sequence
    for enzyme in enz_dict:
        cut_list = []
        count = 0
        ## retrieve cleavage sequence from enzyme dictionary
        cut = enz_dict[enzyme]
        ## find cleavage sites if present in sequence
        p = re.compile(r'(' + cut + ')')
        ## count number of cleavage sites
        it = p.finditer(seq)
        for match in it:
            count += 1
            ## if cleavage sites are found, add to dictionary of cut sequences
            if count != 0:
                cut_list.append((match.start(), match.end()))
                cut_dict[enzyme] = (count, cut_list)

    return cut_dict

#**********************************************************************************

def getEnzyme(acc, enzyme=None):
    """ Function for returning restriction enzyme cleavage sites and indicating 'Bad' or 'Good'.
     Input                      acc                         Gene accession number
                                enzyme                      Optional input for custom cleavage site
     Output                     enzyme_list                 List of all enzymes cutting sequence, Bad/Good,
                                                            and cleavage start/end coordinates
     """

    ## use seq_module to obtain sequence information

    code_seq    = codingSeq(acc)

    base_list = ['A', 'C', 'T', 'G']
    ## call function to show cleavage positions
    if enzyme != None:
        enzyme = enzyme.upper()
        for x in enzyme:
            if x in base_list:
                coding_cut      = enz_cut(acc, code_seq, enzyme)
                seq_cut         = enz_cut(acc, None, enzyme)
            else:
                print('Cleavage site must include A, C, T or G only.')
                break
    else:
        coding_cut      = enz_cut(acc, code_seq)
        seq_cut         = enz_cut(acc)

    ## determine whether enzyme cuts in coding region
    enzymes    = []

    for k,v in seq_cut.items():
        if k in coding_cut:
            enzyme = (k,'Bad', v)
            enzymes.append(enzyme)
        else:
            enzyme = (k,'Good', v)
            enzymes.append(enzyme)
        enzymes.sort()

    return enzymes

#**********************************************************************************


###### main ####

if __name__ == "__main__":


## if file is list of strings, must format sequence:
    # for s in file:
    #     file = s.replace(' ', '')
    #     file = s.upper()
## dummy data

    gene = 'U16860.1'

## get genomic sequence
    line_seq = getSequence(gene)

## print divided into numbered rows of 60
    # for i in range(1, len(line_seq), 60):
    #     print('\n')
    #     print(i, '    ', end='')
    #     for j in range(i, i + 59):
    #         if j in line_seq:
    #             print(line_seq[j], end='')
    # print('\n')


## get annotated sequence with exon boundaries indicated
    with open ('ann_seq_out.txt', 'w') as text_file:

        ann_seq = annotateSeq(gene)
        print(ann_seq, file=text_file)

    text_file.close()
## get coding sequence
    code_seq = codingSeq(gene)
    print(code_seq)

## get genomic sequence with restriction enzyme sites indicated
    custom = 'acnmg'
    enzymes = getEnzyme(gene, custom)
    for x in enzymes:
        print(x)

## get translation
    aa_seq = translate(gene)
    #print(aa_seq[0])
    #print(aa_seq[1])



## line up codons and amino acids
    for x in aa_seq[0]:
        print(x, end=' ')
    print('')
    for x in aa_seq[1]:
        print(x, end='   ')
    print('')




