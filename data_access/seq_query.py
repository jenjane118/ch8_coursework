#!/usr/bin python3

""" Data Access program for sequence information """
"""
Program:        seq_query
File:           seq_query.py

Version:    1.0
Date:       07.05.18
Function:   Query database for sequence for a specific gene using accession number

Author:     Jennifer J. Stiens

Course:     MSc Bioinformatics, Birkbeck University of London
            Biocomputing2 Coursework Assignment

_____________________________________________________________________________

Description:
============
Pymysql query program used to return sequence information for specified gene using accession number.


Usage:
======

seq_query       ACC

Revision History:
=================

v1.0                      07.05.18          Original                    By:Jennifer J. Stiens
                                          
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

#*****************************************************************************

def seq_query(acc):
    """ Return single sequence entry for specified gene.
        Input           acc             accession number
        Output          sequence        (accession number, sequence)
        """

    with cnx.cursor() as cursor:
        # Read a single record
        query = "SELECT accession, sequence FROM sequence WHERE accession = %s;"
        cursor.execute(query, (acc, ))
        sequence = cursor.fetchone()

    return sequence

#*****************************************************************************
## main

if __name__ == "__main__":

    gene = 'AB000381.1'

    seq = seq_query(gene)

    print(seq)