#!/usr/bin python3

""" Data Access program for gene list """

"""
Program:        list_query
File:           list_query.py

Version:    1.0
Date:       07.05.18
Function:   Query database for genebank identifiers for all genes on human chromosome 8

Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/ch8_coursework

Course:     MSc Bioinformatics, Birkbeck University of London
            Biocomputing2 Coursework Assignment

_____________________________________________________________________________

Description:
============
This program will access database using pymysql to obtain list of all genes and gene identifiers from chromosome 8.
Returns tuples containing four identifiers (accession number, gene name, product, chromosomal location).

Usage:
======
list_query


Revision History:
=================

v1.1                     07.05.18           Original         By:Jennifer J. Stiens
                                            
"""

#*****************************************************************************
# Import libraries
import pymysql
from data_access import config_db

#*****************************************************************************

# Connect to MySQL Database
cnx = pymysql.connect(host=config_db.database_config['dbhost'],
                        port=config_db.database_config['port'],
                        user=config_db.database_config['dbuser'],
                        passwd=config_db.database_config['dbpass'],
                        db=config_db.database_config['dbname'])
cursor = cnx.cursor(pymysql.cursors.DictCursor)

def genbank_query():
    with cnx.cursor() as  cursor:
        query = "SELECT accession, gene, product, location FROM genbank;"
        cursor.execute(query)
        genbank = cursor.fetchall()
    return genbank

#*****************************************************************************

###     main    ###

if __name__ == "__main__":
    genes = genbank_query()
    for x in genes:
        print(x)
