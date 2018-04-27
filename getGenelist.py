#!/usr/bin python3

""" Gene Identifiers List"""

"""
Program:    getGenelist
File:       getGenelist.py

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


Usage:
======

Revision History:
=================

V1.0        15.03.18    Original (as Gene_list_class.py)    By: JJS
V1.1        18.03.18    Revised to incorporate methods          JJS
V1.2        22.03.18    included dicttoxml, changed name        JJS

"""

#******************************************************************************
# Import libraries

import gene_module
import pickle, shelve


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
## create gene object for each one
#    gene_object = gene_module.Gene(acc, genid, product, location)
## add to gene_dict using geneList function
#     all_genes = geneList(gene_object)
#
 ## Uses initialisation functi0n to create a gene object for each listing from database


## dummy data here

identity =[('AB061209', 'MRPS28', 'mitochondrial ribosomal protein s28', '8q21.1-q21.2'),
            ('AB098765', 'genid1', 'product1', 'location1'),
            ('AC012345', 'genid2', 'product2', 'location2'),
            ('AC034567', 'genid3', 'product3', 'location3'),
            ('AR456789', 'genid4', 'product4', 'location4')]
gene_list = []
object_dict = {}
## Uses initialisation functi0n to create a gene object for each listing from database
for x in identity:
    gene_object = gene_module.Gene(x[0], x[1], x[2], x[3])
## Adds to list of gene objects
    gene_list.append(gene_object)

## Create a dictionary of gene objects
for object in gene_module.Gene._registry:
    object_dict[object.acc] = (object.genid, object.product, object.location)
for k,v in object_dict.items():
    print(k,v)




xmlTemplate = """
<gene>
<geneidentifiers>
    <accession>%(acc)s</accession>
    <genid>%(genid)s</genid>
    <product>%(product)s</product>
    <location>%(location)s</location>
</geneidentifiers>
</gene>"""



with open('genelist_output.xml', 'w') as xml_file:
    print('<?xml version="1.0" encoding="UTF-8"?>\n', file=xml_file)
    gene_map = {}
    ## Uses gene dictionary function to create a mapping of attributes for printing in xml
    for y in gene_list:
        gene_map = y.geneDict()
        print(xmlTemplate % gene_map, file=xml_file)

## Prints total number of gene objects using static function
print(gene_module.Gene.total())
xml_file.close()



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




