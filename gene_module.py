#!/usr/bin python3

""" Gene Module """

"""
Program:    gene_module
File:       gene_module.py

Version:    1.0
Date:       15.03.18
Function:   Create gene object for genbank database

Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/ch8_coursework
            

Course:     MSc Bioinformatics, Birkbeck University of London
            Biocomputing2 Coursework Assignment

_____________________________________________________________________________

Description:
============
This program includes functions for class of objects called genes

Usage:
======
gene_module         SELF

Revision History:
=================

18.03.18    Revised to incorporate methods      JJS

"""

# ******************************************************************************
# Import libraries


# *****************************************************************************

class Gene:
    """ A simple class for handling associated gene identifiers"""
    count = 0
    _registry = []

    @staticmethod
    def total():
        """ Static method to return current number of gene objects in file"""
        return Gene.count

    # *************************************************************************************

    def __init__(self, acc=str(), genid=str(), product=str(), location=str()):
        """ Create gene object and assigns attributes (identifiers)"""

        self.number = Gene.count
        self.acc = acc
        self.genid = genid
        self.product = product
        self.location = location
        self._registry.append(self)

        Gene.count += 1


    # **************************************************************************************

    def __str__(self):
        """ Print gene object and identifiers"""

        rep = ''
        rep += 'acc: ' + self.acc + '\n' + 'genid: ' + self.genid + '\n' + 'product: ' \
               + self.product + '\n' + 'location: ' + self.location
        return rep

    # **************************************************************************************

    def geneList(self):
        """Create dictionary of gene identifiers"""

        gene_dict = {self.acc: (self.genid, self.product, self.location)}
        return gene_dict

    # **************************************************************************************



#if __name__== "__main__":

