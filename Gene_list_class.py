#!/usr/bin python3

""" Gene Identifiers List"""

"""
Program:    Gene_list_class
File:       Gene_list_class.py

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

#******************************************************************************
# Import libraries

import os
import GeneModule
import pickle, shelve
import dicttoxml

#*****************************************************************************

## dummy data here--in reality would need to import list of gene identifiers
## ideally in lists of tuples for indexing

identity =[('AB061209', 'MRPS28', 'mitochondrial ribosomal protein s28', '8q21.1-q21.2'),
            ('AB098765', 'genid1', 'product1', 'location1'),
            ('AC012345', 'genid2', 'product2', 'location2'),
            ('AC034567', 'genid3', 'product3', 'location3'),
            ('AR456789', 'genid4', 'product4', 'location4')]

#*************************************************************************************

def GeneList(self):
    """Creates a dictionary of gene identifiers.
    """

    gene_dict = {self.acc: (self.genid, self.product, self.location)}
    return gene_dict

#***********************************************************************************
# Main Program
xmlTemplate = """
<gene>
<geneidentifiers>
    <accession>%(acc)s</accession>
    <genid>%(genid)s</genid>
    <product>%(product)s</product>
    <location>%(location)s</location>
</geneidentifiers>
</gene>"""


gene_list = []
all_genes = {}
genes = {}

with open('genelist_output.xml', 'w') as xml_file:
    print('<?xml version="1.0" encoding="UTF-8"?>\n', file=xml_file)

    gene_list = []
    gene_map = {}

    ## Uses initialisation functin to create a gene object for each listing from database
    ## Had to put into a list as function won't work straight from tuple
    for x in identity:
        gene_list.append(GeneModule.Gene(x[0], x[1], x[2], x[3]))

    ## Uses gene dictionary function to create a mapping of attributes for printing in xml
    for y in gene_list:
        gene_map = y.Gene_dict()
        print(xmlTemplate % gene_map, file=xml_file)

## Prints total number of gene objects using static function
print(GeneModule.Gene.total())
xml_file.close()

#for k in genes:
     #   all_genes[k] = genes[k]
#print(all_genes)


#genes_xml = dicttoxml.dicttoxml(all_genes, attr_type=False, custom_root='gene_list')
#xml_file = open('genelist_output.xml', 'w')
#print(genes_xml, file=xml_file)
#xml_file.close()


##****************************************************************
## Other processes:

## Prints total number of gene objects using static function
#print(Gene.total())
## Print gene identifiers using string method
#for y in gene_list:
#    print(y)

## Writes pickled data (dictionary) to .dat file
#f = open('genelist_output.dat', 'wb')
#pickle.dump(gene_list, f)
#f.close()




