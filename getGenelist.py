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
"""

#******************************************************************************
# Import libraries

import gene_module
import sys
from xml.dom import minidom

# Main Program


##  once we have database parser up and running will import python object (list or dictionary)
# from pymysql program import genbank (dictionaries returned: accession:lalala, genid:lalal, etc?)
#all_genes = {}
#for gene in genbank:
#    acc         = gene['accession']
#    genid       = gene['genid']
#    product     = gene['product']
#    location    = gene['location']
#
## Uses initialisation functi0n to create a gene object for each listing from database
#    gene_object = gene_module.Gene(acc, genid, product, location)
#
## dummy data here

identity =[('AB061209', 'MRPS28', 'mitochondrial ribosomal protein s28', '8q21.1-q21.2'),
            ('AB098765', 'genid1', 'product1', 'location1'),
            ('AC012345', 'genid2', 'product2', 'location2'),
            ('AC034567', 'genid3', 'product3', 'location3'),
            ('AR456789', 'genid4', 'product4', 'location4')]

object_dict = {}        #individual object dictionary of identifiers
chrom_dict = {}         #chromosome dictionary of all gene objects

## Uses initialisation function to create a gene object for each listing from database
for x in identity:
    gene_object = gene_module.Gene(x[0], x[1], x[2], x[3])

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

#doc.writexml(sys.stdout, addindent='    ', newl='\n')

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




