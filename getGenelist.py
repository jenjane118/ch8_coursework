#!/usr/bin python3

""" Gene Identifiers List"""

"""
Program:    getGenelist
File:       getGenelist.py

Version:    1.0
Date:       15.03.18
Function:   Retrieve all gene identifier info for every gene from database and returns as a list

Author:     Jennifer J. Stiens

Course:     MSc Bioinformatics, Birkbeck University of London
            Biocomputing2 Coursework Assignment

_____________________________________________________________________________

Description:
============
This program returns a list of all gene identifiers associated with every gene located
on Human Chromosome 8 for display on website.


Usage:
======

Revision History:
=================

V1.0        15.03.18    Original (as Gene_list_class.py)    By: JJS
V1.1        18.03.18    Revised to incorporate methods          JJS
V1.2        22.03.18    included dicttoxml, changed name        JJS
V1.3        29.04.18    changed xml format                      JJS
V1.4        07.05.18    import data access file                 JJS
"""

#******************************************************************************
# Import libraries

import gene_module
import sys
from xml.dom import minidom
from data_access import list_query
#from data_access import config_db

#******************************************************************************
## Main Program ##

genbank = list_query.genbank_query()

for gene in genbank:
    acc         = gene[0]
    genid       = gene[1]
    product     = gene[2]
    location    = gene[3]

## Uses initialisation function to create a gene object for each listing from database
    gene_object = gene_module.Gene(acc, genid, product, location)


object_dict = {}         #individual object dictionary of identifiers
chrom_dict  = {}         #chromosome dictionary of all gene objects

## Create a dictionary of gene objects
for object in gene_module.Gene._registry:
    object_dict = gene_module.Gene.geneList(object)
    for k,v in object_dict.items():
        chrom_dict[k] = v

## write output in xml to file
doc = minidom.Document()

gene = doc.createElement('gene')
doc.appendChild(gene)

for k,v in chrom_dict.items():

    gene_identifiers = doc.createElement('gene_identifiers')
    gene.appendChild(gene_identifiers)

    accession = doc.createElement('accession')
    text = doc.createTextNode(k)
    accession.appendChild(text)
    gene_identifiers.appendChild(accession)

    gene_id = doc.createElement('gene_id')
    text = doc.createTextNode(v[0])
    gene_id.appendChild(text)
    gene_identifiers.appendChild(gene_id)

    product = doc.createElement('product')
    text = doc.createTextNode(v[1])
    product.appendChild(text)
    gene_identifiers.appendChild(product)

    location = doc.createElement('location')
    text = doc.createTextNode(v[2])
    location.appendChild(text)
    gene_identifiers.appendChild(location)

doc.writexml(sys.stdout, addindent='    ', newl='\n')

file_handle = open('genelist_out.xml', 'w')
doc.writexml(file_handle, addindent='   ',newl='\n')
file_handle.close()




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




