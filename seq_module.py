#!/usr/bin python3

""" Sequence annotations module """

"""
Program:        Sequence_module
File:           Sequence_module.py

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
V1.0            Original                        By: JJS
V1.1            Revised as module                   JJS

"""
#*****************************************************************************
# Import libraries

import re



#****************************************************************************


def getSequence(acc, seq):
    """ Returns genomic DNA sequence in numbered form.
    Input               acc                 Accession ID
                        seq                 Sequence

    Output              num_seq             Tuples of line number and 60 bp of seq
    """

    for s in seq:
        dna = s.replace(' ', '')
    num_seq = {}
    count = 0
    for x in dna:
        count += 1
        num_seq[count] = x
    return num_seq

def annotateSeq(acc, seq, exon_list=None):
    """Returns sequence with exon boundaries marked out with symbols
    Input               acc                 Accession ID
                        seq                 Sequence
                        exon_list           List of exon boundaries

    Output              seq                 Annotated sequence with inserted <exon> boundaries
    """

    for s in seq:
        exon_seq = s.replace(' ','')
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

    for s in seq:
        seq = s.replace(' ', '')
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


