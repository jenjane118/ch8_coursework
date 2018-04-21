#!/usr/bin python3

""" Gene Identifiers List"""

"""
Program:    GeneModule
File:       GeneModule.py

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
This program returns a list of all gene identifiers associated with every gene located
on Human Chromosome 8 for display on website.

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
    names = []

    @staticmethod
    def total():
        """ Static method to return current number of gene objects in file.
        """
        return Gene.count

    # *************************************************************************************

    def __init__(self, acc=str(), genid=str(), product=str(), location=str()):
        """ Identifies the gene object and its identifiers
        """
        self.number = Gene.count
        self.acc = acc
        self.genid = genid
        self.product = product
        self.location = location
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

    def Gene_dict(self):
        """Creates a dictionary for mapping to xml template.
        """

        gene_dict = {'acc': self.acc, 'genid': self.genid, 'product': self.product, 'location': self.location}
        return gene_dict

# **************************************************************************************
