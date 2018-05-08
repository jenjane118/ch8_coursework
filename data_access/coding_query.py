#!/usr/bin python3

""" Data Access program for coding region information """
"""
Program:        coding_query
File:           coding_query.py

Version:    1.0
Date:       07.05.18
Function:   Query database for coding information (codon start and exon boundaries) for specific gene using accession number

Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/ch8_coursework

Course:     MSc Bioinformatics, Birkbeck University of London
            Biocomputing2 Coursework Assignment

_____________________________________________________________________________

Description:
============
Pymysql query program used to return coding information for specified gene using accession number.


Usage:
======

coding_query       ACC

Revision History:
=================

v1.0                      07.05.18          Original                                        By:Jennifer J. Stiens

"""
# *****************************************************************************
# Import libraries

import pymysql
from data_access import config_db

# *****************************************************************************

# Connect to MySQL Database
cnx = pymysql.connect(host=config_db.database_config['dbhost'],
                      port=config_db.database_config['port'],
                      user=config_db.database_config['dbuser'],
                      passwd=config_db.database_config['dbpass'],
                      db=config_db.database_config['dbname'])
cursor = cnx.cursor(pymysql.cursors.DictCursor)


# *****************************************************************************

def coding_query(acc):
    """ Return single sequence entry for specified gene.
        Input           acc                 accession number
        Output          coding_info         (accession number, codon start, exon boundaries)
        """

    with cnx.cursor() as cursor:
        # Read a single record
        query = "SELECT accession,  codon_start, positions FROM coding_regions WHERE accession = %s;"
        cursor.execute(query, (acc, ))
        coding_regions = cursor.fetchone()

    return coding_regions


# *****************************************************************************
## main

if __name__ == "__main__":
    gene = 'AB000381.1'

    exons = coding_query(gene)

    print(exons)