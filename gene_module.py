#!/usr/bin python3

""" Gene Identifiers List"""

"""
Program:    gene_module
File:       gene_module.py

Version:    1.0
Date:       15.03.18
Function:   Retrieve all gene identifier info for every gene on
            chromosome 8 and returns as a list

Author:     Jennifer J. Stiens

Course:     MSc Bioinformatics, Birkbeck University of London
            Biocomputing2 Coursework Assignment

_____________________________________________________________________________

Description:
============
This program includes all functions for class of objects called genes

Eventually, will link to a call for detail page for individual genes

Usage:
======

Revision History:
=================

18.03.18    Revised to incorporate methods      JJS

"""

# ******************************************************************************
# Import libraries


# *****************************************************************************


xmlTemplate = """
<gene>
<geneidentifiers>
    <accession>%(acc)s</accession>
    <genid>%(genid)s</genid>
    <product>%(product)s</product>
    <location>%(location)s</location>
</geneidentifiers>
</gene>"""

## dummy data here--in reality would need to import list of gene identifiers
## ideally in lists of tuples for indexing

identity = [('AB061209', 'MRPS28', 'mitochondrial ribosomal protein s28', '8q21.1-q21.2'),
            ('AB098765', 'genid1', 'product1', 'location1'),
            ('AC012345', 'genid2', 'product2', 'location2')]


# *************************************************************************************
class Gene:
    """ A class for handling associated gene identifiers and sequence data."""
    count = 0
    _registry = []

    @staticmethod
    def total():
        """ Static method to return current number of gene objects in file.
        """
        return Gene.count

    # *************************************************************************************

    def __init__(self, acc=str(), genid=str(), product=str(), location=str()):
        """ Creates the gene object and assigns attributes (identifiers)
        """
        self.number = Gene.count
        self.acc = acc
        self.genid = genid
        self.product = product
        self.location = location
        self._registry.append(self)

        Gene.count += 1


    # **************************************************************************************

    def __str__(self):
        """ Function for printing the gene object and identifiers.
        """
        rep = ''
        rep += 'acc: ' + self.acc + '\n' + 'genid: ' + self.genid + '\n' + 'product: ' \
               + self.product + '\n' + 'location: ' + self.location
        return rep

    # ***************************************************************************************

    def geneDict(self):
        """Creates a dictionary for mapping to xml template.
        """

        gene_dict = {'acc': self.acc, 'genid': self.genid, 'product': self.product, 'location': self.location}
        return gene_dict

    # **************************************************************************************

    def geneList(self):
        """Creates a dictionary of gene identifiers.
        """

        gene_dict = {self.acc: (self.genid, self.product, self.location)}
        return gene_dict

    # **************************************************************************************

#if __name__== "__main__":