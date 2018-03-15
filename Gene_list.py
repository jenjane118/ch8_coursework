#!/usr/bin python3

""" Gene Identifiers List"""

"""
Program:    Gene_list
File:       Gene_list.py

Version:    dummy
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



"""

#******************************************************************************
# Import libraries




#*****************************************************************************

xmlTemplate = """<gene>
<geneidentifiers>
    <accession>%(acc)s</accession>
    <genid>%(genid)s</genid>
    <product>%(product)s</product>
    <location>%(location)s</location>
</geneidentifiers>
</gene>"""

## dummy data here--in reality would need to import list of gene identifiers
## ideally in tuples for indexing, struggling to iterate through dictionary
##      and select from multiple values

identity =[('AB061209', 'MRPS28', 'mitochondrial ribosomal protein s28', '8q21.1-q21.2'),
            ('AB098765', 'genid1', 'product1', 'location1'),
            ('AC012345', 'genid2', 'product2', 'location2')]


with open("genelist_output.txt", "w") as text_file:
    for x in identity:   
        gene_list={'acc':x[0], 'genid':x[1], 'product':x[2], 'location':x[3]}
        print(xmlTemplate%gene_list, file=text_file)

 


