#!/usr/bin python3


""" Gene translation Program """

"""
Program:        getEnzymes
File:           getEnzymes.py

Version:        1.0
Date:           21.04.18
Function:       Find and display restriction enzyme cleavage sites in/out of coding region.

Author:         Jennifer J Stiens

Course:         MSc Bioinformatics, Birkbeck University of London
                Biocomputing 2 Coursework

______________________________________________________________________________

Description:
============
This program returns sequence with coding regions and specific restriction enzyme cleavage sites indicated.
User can specify a particular restriction enzyme site or use a custom site. Output is text file in xml format.


Usage:
======

Revision History:
=================
V1.0           21.04.18             Original            By:JJS

"""
#*****************************************************************************
# Import libraries

import sys
import seq_module
from xml.dom import minidom

#****************************************************************************


# dummy data for testing
acc  = 'AB12345'
seq1 = "ATGGCAATGCAGTGGCGCTGTGTCGGACCCGTGCTGTGGCTGCCGAGAGCCATTTTCTGCGAGTGTTTCTCTTCTTCAGGCCCTTTCGGGGTGTAGGCACTGAGAGTGGATCCGAAAGTGGTAGTTCCAATGCCAAGGAGCCTAAGACGCGCGCAGGCGGTTTCGCGAGCGCGTTGGAGCGGCACTCGGAGCTTCTACAGAAGGTGGAGCCCCTACAGAAGGGTTCTCCAAAAAATGTGGAATCCTTTGCATCTATGCTGAGACATTCTCCTCTTACACAGATGGGACCTGCAAAGGATAAACTGGTCATTGGACGGATCTTTCATATTGTGGAGAATGATCTGTACATAGATTTTGGTGGAAAGTTTCATTGTGTATGTAGAAGACCAGAAGTGGATGGAGAGAAATACCAGAAAGGAACCAGGGTCCGGTTGCGGCTATTAGATCTTGAACTTACGTCTAGGTTCCTGGGAGCAACAACAGATACAACTGTACTAGAGGCTAATGCAGTTCTCTTGGGAATCCAGGAGAGTAAAGACTCAAGATCGAAAGAAGAACATCATGAAAAAT"
exons = [('AB12345', 55, 257)]
enz = 'GTGGCG'



## call function to generate exon-annotated sequence
exon_seq    = seq_module.annotateSeq(acc, seq1, exons)

## call function to show cleavage sites in respect to exon boundaries (enz = enzyme name or custom cleavage site)
exon_cut    = seq_module.enz_cut(acc, exon_seq, enz)

## write output in xml to file
doc = minidom.Document()
seq_data = doc.createElement('seq_data')
doc.appendChild(seq_data)

for k,v in exon_cut.items():
    enz = k
    no = str(v[0])
    new_seq = v[1]

    gene = doc.createElement('gene')
    gene.setAttribute('acc', acc)

    seq_data.appendChild(gene)

    sequence = doc.createElement('sequence')
    gene.appendChild(sequence)

    enzyme = doc.createElement('enzyme')
    enzyme.setAttribute('name', enz)
    enzyme.setAttribute('number', no)

    sequence.appendChild(enzyme)

    cut_site = doc.createElement('cut_seq')
    text = doc.createTextNode(new_seq)
    cut_site.appendChild(text)
    enzyme.appendChild(cut_site)

doc.writexml(sys.stdout, addindent='    ', newl='\n')

file_handle = open('getEnzymes_out.xml', 'w')
doc.writexml(file_handle, addindent='   ',newl='\n')
file_handle.close()


