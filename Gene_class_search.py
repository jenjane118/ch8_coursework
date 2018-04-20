#!/usr/bin python3

""" Gene Search Tool"""

"""
Program:    Gene_class_Search
File:       Gene_class_Search.py

Version:    1.0
Date:       18.03.18
Function:   Search list of chromosome 8 genes for one of the identifiers
            and return the other identifiers (and link to details page)

Author:     Jennifer J. Stiens

Course:     MSc Bioninformatics, Birkbeck University of London
            Biocomputing2 Coursework Assignment

_____________________________________________________________________________

Description:
============
This program returns a list of all gene identifiers when searched by acc,
genid, protein product or chromosome location.

Eventually, will link to a call for detail page for individual genes

Usage:
======

Revision History:
=================
V1.0            18.03.18                Original                    By: JJS
V1.1            18.04.18                added 'main'                    JJS


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
 #   list_of_genes = []
    
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
 #       Gene.list_of_genes.append(self)
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

    def Search_acc(list_of_genes, acc):
        """ Search list of genes by accession code and return the
            other identifiers (genid, product, location).

            Input:      list_of_genes           list of all gene objects created
                        acc                     accession number
            Output:     genid                   Gene id
                        product                 protein product
                        location                chromosomal location
        """
        for gene_instance in list_of_genes:
            if acc in gene_instance.acc:
                return(gene_instance)
            
    def Search_genid(list_of_genes, genid):
        """Search list of genes by the gene id and return
            accession number, protein product and location.

            Input:      list_of_genes           list of all gene objects created
                        genid                   Gene id
            Output:     acc                     Accession number
                        genid                   Gene id
                        product                 protein product
                        location                chromosomal location
        """
        for gene_instance in list_of_genes:
            if genid in gene_instance.genid:
                return(gene_instance.acc, gene_instance.genid, gene_instance.product, gene_instance.location)

    def Search_product(list_of_genes, product):
        """Search list of genes by the protein product and return
            accession number, protein product and location.

            Input:      list_of_genes           list of all gene objects created
                        product                 protein product
            Output:     acc                     Accession number
                        genid                   Gene id
                        product                 protein product
                        location                chromosomal location
        """
        for gene_instance in list_of_genes:
            if product in gene_instance.product:
                return(gene_instance.acc, gene_instance.genid, gene_instance.product, gene_instance.location)

    def Search_location(list_of_genes, location):
        """Search list of genes by the chromosomal location and returns
            accession number, gene id, protein product and location.

            Input:      list_of_genes           list of all gene objects created
                        location                chromosomal location
            Output:     acc                     Accession number
                        genid                   Gene id
                        product                 protein product
                        location                chromosomal location
        """
        for gene_instance in list_of_genes:
            if location in gene_instance.location:
                return(gene_instance.acc, gene_instance.genid, gene_instance.product, gene_instance.location)

# main
if __name__ == "__main__":
## Writes output to file in xml format
    with open('genesearch_output.txt', 'w') as text_file:
        print('<?xml version="1.0" encoding="UTF-8"?>\n', file=text_file)
    
        gene_list = []
        search_results = []
        search_map = {}

## Uses initialisation function to create a gene object for each listing from database
## Create a list of gene objects    
        for x in identity:
            gene_list.append(Gene(x[0], x[1], x[2], x[3]))
        
## Uses gene dictionary function to create a mapping of attributes for printing in xml
        search_results = Gene.Search_acc(gene_list, 'AB098765')
        search_map = search_results.Gene_dict()
        print(xmlTemplate%search_map, file=text_file)
    text_file.close()








        






## alternate script without using gene_dict function
## indexes tuple directly (simpler?)
## but then are the attributes assigned?
##with open('genelist_output.txt', 'w') as text_file:
##    for x in identity:
##        gene = Gene(x)
##        gene_map = {'acc': x[0], 'genid': x[1], 'product': x[2], 'location': x[3]}
##        print(xmlTemplate%gene_map, file=text_file)


#### here i directly accessed object's attributes outside of its class definition
##   --not recommended
##        gene_map = {'acc': gene_list[j].acc, 'genid': gene_list[j].genid, 'product': gene_list[j].product, 'location': gene_list[j].location}
##        print(xmlTemplate%gene_map, file=text_file)
##
##
