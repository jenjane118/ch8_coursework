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
This program returns number and positin of restriction enzyme cleavage sites, and indicates whether or not they
lie within the coding region ('usable' or 'bad' tags).
Will test 5 commonly used restriction sites or user can input a custom cleavage site. 
Output is text file in xml format.


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
seq1 = 'agctctcttttttttcccccccgatcgctagctcgctttcgcgcgagaaaaagggctcgcgctagagctcgcgcgggatcggccgtagataattttctctcccccccgcggc'
exons = [('AB12345', 15, 100), ('AB12345', 130, 200)]
enz = 'AAGACCAGAAG'



with open('seq_file2.txt', 'r') as f:
    file = f.read().splitlines()
f.close()

genomic     = seq_module.annotateSeq(acc, file, exons)
print(genomic)
code_seq    = seq_module.codingSeq(acc, file, exons)
print(code_seq)

## call function to show cleavage positions
coding_cut      = seq_module.enz_cut(acc, code_seq, enz)
seq_cut         = seq_module.enz_cut(acc, genomic, enz)
## determine whether enzyme cuts in coding region
good_enzymes    = {}
for k,v in seq_cut.items():
    if k in coding_cut:
        bad_enzymes = {k:v}
        print(k,v)
    else:
        usable = 'Y'
        good_enzymes = {k:v}
        print(k,v)


## write output in xml to file
doc = minidom.Document()
seq_data = doc.createElement('seq_data')
doc.appendChild(seq_data)

for k,v in good_enzymes.items():
    enz = k
    no = str(v[0])
    position = str(v[1])

    gene = doc.createElement('gene')
    gene.setAttribute('acc', acc)

    seq_data.appendChild(gene)

    usable = doc.createElement('usable')
    gene.appendChild(usable)

    enzyme = doc.createElement('enzyme')
    enzyme.setAttribute('name', enz)
    enzyme.setAttribute('cuts', no)

    usable.appendChild(enzyme)

    cut_site = doc.createElement('position')
    text = doc.createTextNode(position)
    cut_site.appendChild(text)
    enzyme.appendChild(cut_site)

for k,v in bad_enzymes.items():
    gene = doc.createElement('gene')
    gene.setAttribute('acc', acc)

    seq_data.appendChild(gene)

    bad = doc.createElement('bad')
    gene.appendChild(bad)

    enzyme = doc.createElement('enzyme')
    enzyme.setAttribute('name', enz)
    enzyme.setAttribute('cuts', no)

    bad.appendChild(enzyme)

    cut_site = doc.createElement('position')
    text = doc.createTextNode(position)
    cut_site.appendChild(text)
    enzyme.appendChild(cut_site)

doc.writexml(sys.stdout, addindent='    ', newl='\n')

file_handle = open('getEnzymes_out.xml', 'w')
doc.writexml(file_handle, addindent='   ',newl='\n')
file_handle.close()


