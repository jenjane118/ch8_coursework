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
V1.5            07.05.18    Linking to data access scripts      JJS
                            Changing coding info scripts        JJS
"""
#*****************************************************************************
# Import libraries


import re
import sys
from data_access import seq_query
from data_access import coding_query

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

#**********************************************************************************

def getSeq(acc):
    """Query database for  genomic sequence and test for valid sequence.
    Input           acc                 Accession ID
    Output          gen_seq             Genomic dna sequence
    """

    ## use pymysql script with accession number to access database for 'sequence'

    gen_seq = seq_query.seq_query(acc)


    # ## check to see accession number requested matches return
    # ## extract sequence
    try:
        gene = gen_seq[0]
    except TypeError:
        print('Gene not found in database')
        exit(0)
    if gene == acc:
        seq = gen_seq[1]
    else:
        print('Error: Sequence not found')
        exit(0)

    ## check to see if sequence is valid (made up of valid nucleotide symbols (a,c,t,g,n)
    nucleotides = {'a': True, 'c': True, 't': True, 'g': True, 'n': True}
    for letter in seq:
        try:
            nucleotides[letter]
        except KeyError:
            print('Sequence not valid')
            exit(0)

    return gen_seq

#**********************************************************************************

def getCoding(acc):
    """Return genomic sequence, codon start, exon positions.
        Input           acc                 Accession ID
        Output          coding_dict         {acc: codon start, exon positions}

        """

    ## use getCoding(acc) function to get code start and exon boundaries from pymysql script
    code_info   = coding_query.coding_query(acc)
    gene        = code_info[0]
    code_start  = code_info[1]
    positions   = code_info[2]

    ## use regex to extract exon boundaries
    p = re.compile(r'\d+')
    exon_list = p.findall(positions)

    ## make list of start/end pairs
    new_list = []
    for i in range(0, len(exon_list),2):
        exon_pair  = (int(exon_list[i]), int(exon_list[i+1]))
        new_list.append(exon_pair)

    coding_dict = {gene: (code_start, new_list)}

    return coding_dict

#****************************************************************************

def numSequence(acc):
    """ Return genomic DNA sequence in numbered form.
    Input               acc                 Accession ID

    Output              num_seq             Dictionary of numbered bp for requested sequence
    """
    ## use function, getSeq(acc) to get sequence info from pymysql script
    seq_info = getSeq(acc)
    acc = seq_info[0]
    seq = seq_info[1]

    ## format sequence in unbroken, uppercase form
    seq = seq.replace(' ', '')
    seq = seq.upper()

    num_seq = {}
    count = 0
    for x in seq:
        count += 1
        num_seq[count] = x
    return num_seq

#**********************************************************************************

def annotateSeq(acc):
    """Return sequence with exon boundaries marked out with symbols (or simple sequence string if no exons indicated).
    Input               acc                 Accession ID

    Output              exon_seq            Annotated sequence with inserted *exon/exon* boundaries (string)
    """

    ## use function, getSeq(acc) to get sequence info from pymysql script
    seq_info    = getSeq(acc)
    seq_gene    = seq_info[0]
    seq         = seq_info[1]

    ## format sequence into unbroken uppercase string
    seq = seq.replace(' ', '')
    seq = seq.upper()

    ## use getCoding(acc) function to get code start and exon boundaries from pymysql script
    coding = getCoding(acc)
    try:
        coding[seq_gene]
        exon_list = coding[seq_gene][1]
    except KeyError:
        print('Gene not found.')
        exit(0)

    ## create new sequence string and identify index for exon start and end
    exon_seq        = ''
    base_count      = 0
    if exon_list   != None:
        for base in seq:
            base_count  += 1
            exon_seq    += base
            for pair in exon_list:
                exon_start = (pair[0]) -1           #correction for indexing
                exon_end = pair[1]
            ## insert 'exon' markers at start/end of exons
                if base_count == exon_start:
                    exon_seq += '*exon'
                elif base_count == exon_end:
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

    ## use function, getSeq(acc) to get sequence info from pymysql script
    seq_info = getSeq(acc)
    seq_gene = seq_info[0]
    seq      = seq_info[1]

    ## reformat sequence into unbroken uppercase string
    seq = seq.replace(' ', '')
    seq = seq.upper()

    ## use getCoding(acc) function to get code start and exon boundaries from pymysql script
    coding      = getCoding(acc)
    try:
        coding[seq_gene]
        codon_start = coding[seq_gene][0]
        exon_list   = coding[seq_gene][1]
    except KeyError:
        print('Gene not found.')
        exit(0)

    ## assemble coding sequence by using exon boundaries as index start/end
    if exon_list       != None:
        coding_seq      = ''
        first_exon      = True
        # have to subtract one from codon_start as index correction
        codon_start     = codon_start - 1
        ## first exon needs to use codon start to initiate translation
        for x in exon_list:
            if first_exon   == True:
                start       = (x[0]) - 1          # must subtract 1 for index correction
                start       = start + codon_start
                end         = x[1]
                first_exon  = False
            else:
                start   = (x[0]) - 1
                end     = x[1]
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
    codon_list  = []
    codon       = ''

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
        ## use function, getSeq(acc) to get sequence info from pymysql script
        seq_info = getSeq(acc)
        acc = seq_info[0]
        seq = seq_info[1]

    ## reformat sequence into unbroken upperclass string
    seq = seq.replace(' ', '')
    seq = seq.upper()

    cut_dict = {}
    ## if a custom cleavage site is included, will indicate number and position of cleavage sites in sequence
    if enzyme != None:
        cut_list    = []
        count       = 0
        cut         = enzyme

        ## find custom cleavage sites if present in sequence
        p = re.compile(r'(' + cut + ')')
        it = p.finditer(seq)

        ## count number of cleavage sites
        for match in it:
            count += 1
            ## if cleavage site is present in sequence, add number of cuts and site location to cut_dict
            if count != 0:
                cut_list.append((match.start(), match.end()))
                cut_dict[enzyme] = (count, cut_list)

    ## for each enzyme in dictionary, it will search and return only the enzymes which have cleavage sites in sequence
    for enzyme in enz_dict:
        cut_list    = []
        count       = 0

        ## retrieve cleavage sequence from enzyme dictionary
        cut = enz_dict[enzyme]

        ## find cleavage sites if present in sequence
        p = re.compile(r'(' + cut + ')')
        ## count number of cleavage sites
        it = p.finditer(seq)
        for match in it:
            count += 1
            ## if cleavage sites are found, add number of cuts and site location to cut_dict
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
    ## make dictionary for easier searching of enzymes
    results_dict = {}
    for x in enzymes:
        results_dict[x[0]] = (x[1], x[2])


    return results_dict

#**********************************************************************************


###### main ####

if __name__ == "__main__":

    gene = 'AB000381.1'
    coding = getCoding(gene)
    print(coding)

## get genomic sequence
    line_seq = numSequence(gene)

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

## get coding sequence
    code_seq = codingSeq(gene)
    print(code_seq)

## get genomic sequence with restriction enzyme sites indicated
    custom = 'gggg'
    enzymes = getEnzyme(gene, enzyme = 'gggg')
    for k,v in enzymes.items():
        print(k, ':',v)

## get translation
    aa_seq = translate(gene)

## line up codons and amino acids
    for x in aa_seq[0]:
        print(x, end=' ')
    print('')
    for x in aa_seq[1]:
        print(x, end='   ')
    print('')




