#!/usr/bin python3

""" Gene Identifiers List"""

"""
Program:    Gene_list
File:       Gene_list.py

Version:    1.0
Date:       15.03.18
Function:   Retrieve all gene identifier info for every gene on
            chromosome 8 and returns as a list

Author:     Jennifer J. Stiens

Course:     MSc Bioninformatics, Birkbeck University of London
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

#******************************************************************************
# Import libraries

import os

#*****************************************************************************


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

identity =[('AB061209', 'MRPS28', 'mitochondrial ribosomal protein s28', '8q21.1-q21.2'),
            ('AB098765', 'genid1', 'product1', 'location1'),
            ('AC012345', 'genid2', 'product2', 'location2')]


class Gene:
    """ A class for handling associated gene identifiers and sequence data."""
    count = 0
    names = []

    @staticmethod
    def total():
        """ Static method to return current number of gene objects in file.
        """
        return Gene.count
    
    def __init__(self, acc= str() , genid = str(), product=str(), location=str()):
        self.number     = Gene.count
        self.acc        = acc
        self.genid      = genid
        self.product    = product
        self.location   = location
        Gene.count      += 1

    def __str__(self):
        rep = ''
        rep += 'acc: ' + self.acc + '\n' + 'genid: ' + self.genid + '\n' + 'product: ' \
               + self.product  + '\n' + 'location: ' + self.location
        return rep
 
    def Gene_dict(self):
        """Create a dictionary for mapping to xml template.
        """
        
        gene_dict = {'acc': self.acc, 'genid': self.genid, 'product': self.product, 'location': self.location} 
        return gene_dict



## Writes output to file in xml format
with open('genelist_output.txt', 'w') as text_file:
    print('<?xml version="1.0" encoding="UTF-8"?>\n', file=text_file)
    
    gene_list = []
    gene_map = {}

## Uses initialisation functin to create a gene object for each listing from database
## Had to put into a list as function won't work straight from tuple    
    for x in identity:
        gene_list.append(Gene(x[0], x[1], x[2], x[3]))
        
## Uses gene dictionary function to create a mapping of attributes for printing in xml
    for y in gene_list:
        gene_map = y.Gene_dict()
        print(xmlTemplate%gene_map, file=text_file)

## Prints total number of gene objects using static function
print(Gene.total())


