#!/usr/bin python3


""" Restriction Enzyme Cleavage Function """

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
This program returns number and position of restriction enzyme cleavage sites, and indicates whether or not they
lie within the coding region. Returns list of Enzymes, whether they are a 'good' or 'bad' enzyme, number of cleavage
sites, coordinates of cleavage sites.
Will test 5 commonly used restriction enzymes and user can input a custom cleavage site. 


Usage:
======

Revision History:
=================
V1.0           21.04.18             Original                            By:JJS
V1.1           24.04.18             Revised                                JJS
V1.2           26.04.18             Rewritten as function                  JJS 
"""
#*****************************************************************************
# Import libraries

import sys
import seq_module
from xml.dom import minidom

#****************************************************************************

## main ##

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

##  main ##
if __name__ == "__main__":
# dummy data for testing
    acc  = 'AB12345'

    enz = 'GAGAGC'

    cut_sites = getEnzyme(acc, enz)
    print(cut_sites)

