#!/usr/bin python3

""" Sequence annotations module """

"""
Program:        Seq_module
File:           Seq_module.py

Version:        1.0
Date:           22.03.18
Function:       Return genomic sequence information for specified gene.
                Indicate positions of exons.

Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
this program returns genomic sequence information for specified gene.


Usage:
======

Revision History:
=================
V1.0            22.03.18    Original                        By: JJS
V1.1            06.04.18    Revised as module                   JJS
V1.2            21.04.18    Combined translate function         JJS
"""
#*****************************************************************************
# Import libraries


import re


#****************************************************************************


def getSequence(acc, seq):         ## need to enter sequence until get usable data from database
    """ Returns genomic DNA sequence in numbered form.
    Input               acc                 Accession ID
                        seq                 Sequence

    Output              num_seq             Dictionary of numbered bp of sequence
    """
    ## use accession number to retrieve sequence from database
    ## request database to return sequence for specific accession number
    # for x in sequence_file:
    ## check to see accession number requested matches return
    #      if x[0] = acc:
    #          x[1] = seq
    #      else:
    #          print('Error: Accession numbers do not match')

    ## format sequence in unbroken, upperclass form
    seq_upper = ''
    for s in seq:
        seq = s.replace(' ', '')
        seq_upper += s.upper()

    num_seq = {}
    count = 0
    for x in seq_upper:
        count += 1
        num_seq[count] = x
    return num_seq

def annotateSeq(acc, seq, exon_list=None):           ## need to enter sequence, exon_list until get usable data from db
    """Returns sequence with exon boundaries marked out with symbols
    Input               acc                 Accession ID
                        seq                 Sequence
                        exon_list           List of exon boundaries

    Output              exon_seq            Annotated sequence with inserted *exon/exon* boundaries
    """

    ## use accession number to retrieve sequence from database
    ## request database to return sequence for specific accession number
    # for x in exon_file:
    ## check to see accession number requested matches return
    #      if x[0] = acc:
    #
    # x[1] = seq
    #      else:
    #          print('Error: Accession numbers do not match')

    seq_dict = getSequence(acc, seq)

    exon_seq = ''
    # for s in seq:
    #    seq = s.replace(' ','')
    #    exon_seq += s.upper()
    ## in list of info, identify start and end sites for exon boundaries (might change with format from database?)
    if exon_list != None:
        for k in seq_dict:
            for x in exon_list:
                start = x[1]
                end = x[2]
                exon_seq += seq_dict[k]
                if k == start:
                    exon_seq += '*exon'
                elif k == end:
                    exon_seq += 'exon*'
    return exon_seq


def codingSeq(acc, seq, exon_list=None):
    """Returns coding sequence (stuck together exons).
    Input           acc                 Accession ID
                    seq                 Sequence
                    exon_list           List of exon boundaries

    Output          coding_seq          Coding sequence
    """

    #change formatting of string
    for s in seq:
        seq = s.replace(' ', '')
        seq = s.upper()
    if exon_list != None:
        coding_seq = ''
        for x in exon_list:
            start       = x[1]
            end         = x[2]
            #code_start = x[3]
            coding_seq += seq[start:end]
            #coding_seq = coding_seq[code_start:]
        return coding_seq
    else:
        return coding_seq

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

    ## divides string into list of codons
    for x in dna:
        codon += x
        if len(codon) == 3:
            codon_list.append(codon)
            codon = ''

    ## finds corresponding amino acid for codon in list
    for codon in codon_list:
        if codon in codon_table:
            aa_seq += codon_table[codon]
    return codon_list, aa_seq


def enz_cut(acc, sequence, enzyme=None):
    """ Will indicate any cleavage sites from popular restriction enzymes in
        restriction enzyme dictionary. User can also search a custom cleavage site.
        Sequence is returned with indicated cutting sites ('@' symbol) and exon boundaries ('<exon' and 'exon>').

    Input           acc                 accession number for gene
                    sequence            gene sequence (raw)
                    enzyme              optional custom cleavage site

    Output          cut_tuple           tuple containing:
                                        (enzyme/cleavage sequence, no. of cleavage sites, annotated seq)
    """

    enz_dict = {
        'EcoRI': 'GAATTC', 'BamHI': 'GGATCC',
        'BsuMI': 'CTCGAG', 'HindIII': 'AAGCTT',
        'EcoRV': 'GATATC'}

    cut_dict = {}
    ## if a custom cleavage site is included, will indicate number and position of cleavage sites in sequence
    if enzyme != None:
        count = 0
        cut = enzyme
        new_seq = sequence.replace(' ', '')
        ## find custom cleavage sites if present in sequence
        p = re.compile(r'(' + cut + ')')
        it = p.finditer(new_seq)
        ## count number of cleavage sites
        for match in it:
            count += 1
        ## put '@' symbols around the identified cleavage site
        cut_seq = re.sub(r'(' + cut + ')', r'@\1@', new_seq)
        ## if cleavage site is present in sequence, add to cut dictionary
        if count !=0:
            cut_dict[enzyme] = (count, cut_seq)


    ## for each enzyme in dictionary, it will search and return only the enzymes which have cleavage sites in sequence
    for enzyme in enz_dict:
        count = 0
        ## retrieve cleavage sequence from enzyme dictionary
        cut = enz_dict[enzyme]
        new_seq = sequence.replace(' ', '')
        ## find cleavage sites if present in sequence
        p = re.compile(r'(' + cut + ')')
        ## count number of cleavage sites
        it = p.finditer(new_seq)
        for match in it:
            count += 1
        ## place '@' symbols around identified cleavage sites
        cut_seq = re.sub(r'(' + cut + ')', r'@\1@', new_seq)
        ## if cleavage sites are found, add to dictionary of cut sequences
        if count != 0:
            cut_dict[enzyme] = (count, cut_seq)
    return cut_dict

# main
if __name__ == "__main__":

## dummy data (sequence in file, acc and exon_list below)
    with open ('seq_file1.txt', 'r') as f:
        file = f.read().splitlines()
    f.close()

## dummy data
    exon_list = [('AB12345.1', 36, 50), ('AB12345.1', 55, 70)]
    acc = 'AB12345.1'

## get genomic sequence
    line_seq = getSequence(acc, file)
## print without line numbers
    for k, v in line_seq.items():
        print(v, end='')

## print divided into numbered rows of 60
    # for i in range(1, len(line_seq), 60):
    #     print('\n')
    #     print(i, '    ', end='')
    #     for j in range(i, i + 59):
    #         if j in line_seq:
    #             print(line_seq[j], end='')

## get annotated sequence with exon boundaries indicated
    ann_seq = annotateSeq(acc, file, exon_list)
    print('\n', ann_seq)

## get coding sequence
    code_seq = codingSeq(acc, file, exon_list)
    print(code_seq)

## get genomic sequence with restriction enzyme sites (and exons) indicated
    enz_seq = enz_cut(acc, ann_seq)
    for k, v in enz_seq.items():
        print(k,v)

## get translation
    aa_seq = translate(acc, code_seq)

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

    xml_seq_ann = """
    <gene>
        <sequence acc=%(acc)s>
            <ann_seq>%(sequence)s</ann_seq>
        </sequence>
    </gene>
    """


    seq_map = {'acc': acc, 'sequence': ann_seq}

    xml_seq_coding = """
    <gene>
        <sequence acc=%(acc)s>
            <coding_seq>%(coding_seq)s</coding_seq>
        </sequence>
    </gene>
    """
    seq_map2 = {'acc':acc, 'coding_seq': code_seq}

    xml_seq_translation = """
    <gene>
        <sequence acc=%(acc)s>
            <aa_seq>%(aa_seq)s</aa_seq>
        </sequence>
    </gene>
    """
    seq_map3 = {'acc':acc, 'aa_seq':aa_seq[1]}

    xml_seq_codons = """
    <gene>
        <sequence acc=%(acc)s>
            <codon_seq>%(codon_seq)s</codon_seq>
        </sequence>
    </gene>
    """
    seq_map4 = {'acc':acc, 'codon_seq':aa_seq[0]}



    with open('sequence_output.xml', 'w') as text_file:
        print('<?xml version="1.0" encoding="UTF-8"?>', file=text_file)
        print(xml_seq_ann %seq_map, file=text_file)
        print(xml_seq_coding %seq_map2, file=text_file)
        print(xml_seq_translation %seq_map3, file=text_file)
        print(xml_seq_codons %seq_map4, file=text_file)
    text_file.close()


