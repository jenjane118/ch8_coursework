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





#****************************************************************************


def getSequence(acc, seq):
    """ Returns genomic DNA sequence in numbered form.
    Input               acc                 Accession ID
                        seq                 Sequence

    Output              num_seq             Dictionary of numbered bp of sequence
    """
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

def annotateSeq(acc, seq, exon_list=None):
    """Returns sequence with exon boundaries marked out with symbols
    Input               acc                 Accession ID
                        seq                 Sequence
                        exon_list           List of exon boundaries

    Output              exon_seq            Annotated sequence with inserted <exon> boundaries
    """
    exon_seq = ''
    for s in seq:
        seq = s.replace(' ','')
        exon_seq += s.upper()
    if exon_list != None:
        for x in exon_list:
            acc     = x[0]
            start   = x[1]
            end     = x[2]
            exon_seq     = exon_seq[:start] + '<exon' + exon_seq[start:end] +'exon>' + exon_seq[end:]
        return exon_seq


def codingSeq(acc, seq, exon_list=None):
    """Returns coding sequence (stuck together exons).
    Input           acc                 Accession ID
                    seq                 Sequence
                    exon_list           List of exon boundaries

    Output          coding_seq          Coding sequence
    """
    seq_upper = ''
    #change formatting of string
    for s in seq:
        seq         = s.replace(' ', '')
        seq_upper  += s.upper()
    if exon_list != None:
        coding_seq = ''
        for x in exon_list:
            start       = x[1]
            end         = x[2]
            #code_start = x[3]
            coding_seq += seq_upper[start:end]
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

    # divides string into list of codons
    codon_list = []
    codon = ''
    for x in dna:
        codon += x
        if len(codon) == 3:
            codon_list.append(codon)
        for codon in codon_list:            # finds corresponding amino acid for codon in list
            if codon in codon_table:
                aa_seq += codon_table[codon]
        return codon_list, aa_seq


# main
if __name__ == "__main__":

## dummy data (sequence in file, acc and exon_list below)
    with open ('seq_file1.txt', 'r') as f:
        file = f.read().splitlines()
    f.close()

    exon_list = [('AB12345.1', 36, 807), ('AB232010.1', 2222, 2311)]
    acc = 'AB12345.1'

# call functions
    line_seq = getSequence(acc, file)
    #print(line_seq)

    ann_seq = annotateSeq(acc, file, exon_list)
    #print(ann_seq)

    code_seq = codingSeq(acc, file, exon_list)
    #print(code_seq)


# print to file in xml format
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



## # Print out the annotated gene sequence of a particular gene in xml format:
    with open('sequence_output.txt', 'w') as text_file:
    # I can't manage to get numbered sequence to show up in xml--it is in form of dictionary (line_seq).
    # here is script for printing dictionary in tidy rows with line numbers:
        for i in range(1, len(line_seq), 60):
            print(i, '    ', end='', file=text_file)
            for j in range(i, i + 59):
                if j in line_seq:
                    print(line_seq[j], end='', file=text_file)
                else:
                    pass
            print('\n', file=text_file)
        print('<?xml version="1.0" encoding="UTF-8"?>', file=text_file)
        print(xml_seq_ann %seq_map, file=text_file)
        print(xml_seq_coding %seq_map2, file=text_file)
    text_file.close()


