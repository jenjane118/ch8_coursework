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
    seq = 'atggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacgagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaa'

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
    seq = 'atggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacgagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaa'
    exon_list = [('AB12345.1', 36, 150), ('AB12345.1', 155, 270)]

    exon_seq = ''

    seq = seq.replace(' ', '')
    seq = seq.upper()

    count   = 0
    if exon_list    != None:
        for base in seq:
            count += 1
            exon_seq += base
            for x in exon_list:
                start = x[1]-1          ## subtract one for index start
                end = x[2]-1
                ## insert 'exon' markers at start/end of exons
                if count == start:
                    exon_seq += '*exon'
                elif count == end:
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

    ## use pymysqul script with accession number to retrieve 'exon'
    ## check to see accession number requested matches database return
    ## make list of exons
    # for x in exon_file:
    #    if x[0] = acc:
    #        x[1] = exon
    #        exon_list.append(exon)

    ## dummy sequence until have data
    seq = 'atggcggcgctgtgtcggacccgtgctgggggtaggcgcgcgcgcgcgcggtaggcactgagagtggtacatatttttagggatatcgctcgagagagagagacacacacacacaccacaccacaatatattatattatattagggagaggaatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacgagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaa'
    exon_list = [('AB12345.1', 36, 105), ('AB12345.1', 155, 270)]

    seq = seq.replace(' ', '')
    seq = seq.upper()

    if exon_list != None:
        coding_seq = ''
        for x in exon_list:
            start       = x[1]-1        # subtract for index starting at 0
            end         = x[2]-1
            #code_start = x[3]
            coding_seq += seq[start:end]
            #coding_seq = coding_seq[code_start:]
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
    return codon_list, aa_seq

#**********************************************************************************

def enz_cut(acc, seq=None, enzyme=None):
    """ Indicate any cleavage sites from restriction enzymes in
        restriction enzyme dictionary. Optionally, search a custom cleavage site (enzyme).
        Return dictionary of enzyme name and cleavage positions.

    Input           acc                 accession number for gene
                    sequence            gene sequence string (optional)
                    enzyme              custom cleavage site (optional)

    Output          cut_dict            dictionary {enzyme:no. of cleavage sites, positions}

    """

    enz_dict = {
        'EcoRI': 'GAATTC', 'BamHI': 'GGATCC',
        'BsuMI': 'CTCGAG', 'HindIII': 'AAGCTT',
        'EcoRV': 'GATATC'}

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
        seq = 'atggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacgagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaaatggcggcgctgtgtcggacccgtgctgtggctgccgagagccattttctgcgagtgtttctcttcttcaggccctttcggggtgtaggcactgagagtggatccgaaagtggtagttccaatgccaaggagcctaagacgcgcgcaggcggtttcgcgagcgcgttggagcggcactcggagcttctacagaaggtggagcccctacagaagggttctccaaaaaatgtggaatcctttgcatctatgctgagacattctcctcttacacagatgggacctgcaaaggataaactggtcattggacggatctttcatattgtggagaatgatctgtacatagattttggtggaaagtttcattgtgtatgtagaagaccagaagtggatggagagaaataccagaaaggaaccagggtccggttgcggctattagatcttgaacttacgtctaggttcctgggagcaacaacagatacaactgtactagaggctaatgcagttctcttgggaatccaggagagtaaagactcaagatcgaaagaagaacatcatgaaaaa'

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
    genomic     = seq_module.annotateSeq(acc)
    code_seq    = seq_module.codingSeq(acc)


    ## call function to show cleavage positions
    if enzyme != None:
        coding_cut      = seq_module.enz_cut(acc, code_seq, enz)
        seq_cut         = seq_module.enz_cut(acc, genomic, enz)
    else:
        coding_cut      = seq_module.enz_cut(acc, code_seq)
        seq_cut         = seq_module.enz_cut(acc, genomic)

    ## determine whether enzyme cuts in coding region
    enzymes    = []

    for k,v in seq_cut.items():
        if k in coding_cut:
            enzyme = (k,'Bad', v)
            enzymes.append(enzyme)
        else:
            enzyme = (k,'Good', v)
            enzymes.append(enzyme)

    return enzymes

#**********************************************************************************


###### main ####

if __name__ == "__main__":

## dummy data (sequence in file, acc and exon_list below)
    # with open ('AB000381.1.txt', 'r') as f:
    #     file = f.read().splitlines()
    # f.close()

## if file is list of strings, must format sequence:
    # for s in file:
    #     file = s.replace(' ', '')
    #     file = s.upper()
## dummy data
    #exon_list = [('AB12345.1', 36, 50), ('AB12345.1', 55, 70)]
    gene = 'AB000381.1'

## get genomic sequence
    line_seq = getSequence(gene)

## print sequence
    for k, v in sorted(line_seq.items()):
        print(v, end='')
    print('\n')

## print divided into numbered rows of 60
    # for i in range(1, len(line_seq), 60):
    #     print('\n')
    #     print(i, '    ', end='')
    #     for j in range(i, i + 59):
    #         if j in line_seq:
    #             print(line_seq[j], end='')
    # print('\n')

## get annotated sequence with exon boundaries indicated
    ann_seq = annotateSeq(gene)
    print(ann_seq)

## get coding sequence
    code_seq = codingSeq(gene)
    print(code_seq)

## get genomic sequence with restriction enzyme sites indicated (if exons indicated, it can interfere with recognition)

    enz_seq = enz_cut(gene, 'ACTT')
    for k, v in enz_seq.items():
        print(k,v)

## get translation
    aa_seq = translate(gene)
    print(aa_seq[0])
    print(aa_seq[1])



## line up codons and amino acids
#     for x in aa_seq[0]:
#         print(x, end=' ')
#     print('')
#     for x in aa_seq[1]:
#         print(x, end='   ')
#     print('')
# ## print amino acid sequence
#     print(aa_seq[1])


## if you want to print to file in xml format to .xml file:
    #
    # xml_seq_ann = """
    # <gene>
    #     <sequence acc=%(acc)s>
    #         <ann_seq>%(sequence)s</ann_seq>
    #     </sequence>
    # </gene>
    # """
    #
    #
    # seq_map = {'acc': acc, 'sequence': ann_seq}
    #
    # xml_seq_coding = """
    # <gene>
    #     <sequence acc=%(acc)s>
    #         <coding_seq>%(coding_seq)s</coding_seq>
    #     </sequence>
    # </gene>
    # """
    # seq_map2 = {'acc':acc, 'coding_seq': code_seq}
    #
    # xml_seq_translation = """
    # <gene>
    #     <sequence acc=%(acc)s>
    #         <aa_seq>%(aa_seq)s</aa_seq>
    #     </sequence>
    # </gene>
    # """
    # seq_map3 = {'acc':acc, 'aa_seq':aa_seq[1]}
    #
    # xml_seq_codons = """
    # <gene>
    #     <sequence acc=%(acc)s>
    #         <codon_seq>%(codon_seq)s</codon_seq>
    #     </sequence>
    # </gene>
    # """
    # seq_map4 = {'acc':acc, 'codon_seq':aa_seq[0]}
    #


    # with open('sequence_output.xml', 'w') as text_file:
    #     print('<?xml version="1.0" encoding="UTF-8"?>', file=text_file)
    #     print(xml_seq_ann %seq_map, file=text_file)
    #     print(xml_seq_coding %seq_map2, file=text_file)
    #     print(xml_seq_translation %seq_map3, file=text_file)
    #     print(xml_seq_codons %seq_map4, file=text_file)
    # text_file.close()
    #
    #
