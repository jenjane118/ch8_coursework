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

"""
#*****************************************************************************
# Import libraries

import re
import sys
from xml.dom import minidom

#****************************************************************************



def getSequence(acc, seq):
    """ Returns genomic DNA sequence in numbered form.
    Input               acc                 Accession ID
                        seq                 Sequence

    Output              num_seq             Tuples of line number and 60 bp of seq
    """

    for s in seq:
        dna = s.replace(' ', '')
    num_seq = []
    gen_seq = ''
    for i in range(1, len(dna), 60):
        x = (i, dna[i-1:i+59])
        num_seq.append(x)

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
            start   = x[0]
            end     = x[1]
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
        print(len(seq))
    if exon_list != None:
        coding_seq = ''
        for x in exon_list:
            print(x)
            start       = x[0]
            end         = x[1]
            coding_seq += seq[start:end]
            print(coding_seq)
            print(len(coding_seq))
        return coding_seq







xml_seq_ann = """
<gene>
<sequence acc=%(acc)s>
<ann_seq>%(sequence)s</ann_seq>
</sequence>
</gene>
"""

xml_seq_coding = """
<gene>
<sequence acc=%(acc)s>
<coding_seq>%(coding_seq)s</coding_seq>
</sequence>
</gene>
"""


with open ('seq_file1.txt', 'r') as f:
    file = f.read().splitlines()


exon_list = [(36, 807), (2010, 2222)]
acc = 'AB12345'

    
num_seq = getSequence(acc, file)

ann_seq = annotateSeq(acc, file, exon_list)

code_seq = codingSeq(acc, file, exon_list)





seq_map = {'acc': acc, 'sequence': ann_seq}
seq_map2 = {'acc':acc, 'coding_seq': code_seq}



# Print out the annotated gene sequence of a particular gene in xml format:
with open('sequence_output.txt', 'w') as text_file:
    # I can't manage to get numbered sequence to show up in xml--it is in form of tuples. Not sure we need this anyway
    # could return numbered dictionary instead?
    print('\nNumbered genomic sequence:\n', num_seq, file=text_file)

    print('<?xml version="1.0" encoding="UTF-8"?>', file=text_file)
    print(xml_seq_ann %seq_map, file=text_file)
    print(xml_seq_coding %seq_map2, file=text_file)


