#!/usr/bin python3


## dummy file for testing out restriction enzyme cutting sites

import re
from re import finditer



def get_enzyme(enzyme, sequence):
    """Get restriction enzyme cutting sequence
        return cutting sites (base numbers) in dna sequence."""
    cut_list = []
    for match in finditer(enzyme, sequence):
        cut_list.append(match.span())
    return cut_list

    

enz = 'fake enzyme'
fake_cut = 'AATGCAGT'

acc  = 'AB12345'
seq1 = 'ATGGCAATGCAGTGGCGCTGTGTCGGACCCGTGCTGTGGCTGCCGAGAGCCATTTTCTGCGAGTGTTTCTCTTCTTCA\
        GGCCCTTTCGGGGTGTAGGCACTGAGAGTGGATCCGAAAGTGGTAGTTCCAATGCCAAGGAGCCTAAGA\
        CGCGCGCAGGCGGTTTCGCGAGCGCGTTGGAGCGGCACTCGGAGCTTCTACAGAAGGTGGAGCCCCTACA\
        GAAGGGTTCTCCAAAAAATGTGGAATCCTTTGCATCTATGCTGAGACATTCTCCTCTTACACAGATGGGA\
        CCTGCAAAGGATAAACTGGTCATTGGACGGATCTTTCATATTGTGGAGAATGATCTGTACATAGATTTTG\
        GTGGAAAGTTTCATTGTGTATGTAGAAGACCAGAAGTGGATGGAGAGAAATACCAGAAAGGAACCAGGGT\
        CCGGTTGCGGCTATTAGATCTTGAACTTACGTCTAGGTTCCTGGGAGCAACAACAGATACAACTGTACTA\
        GAGGCTAATGCAGTTCTCTTGGGAATCCAGGAGAGTAAAGACTCAAGATCGAAAGAAGAACATCATGAAA\
        AAT'

##  run get_enzyme function and output into list       

span_list = []
span_list = get_enzyme(fake_cut, seq1)

## convert list of match items to strings. Add one to each position (starts at 0). Append to new list
new_list =[]
n_cuts = 0
#print(span_list)
for x in span_list:
    n_cuts += 1
    span   =  str(x[0]+1) + '-' + str(x[1]+1)
    new_list.append(span)
#print(new_list)
#print(n_cuts)


## print output in xml to file

from xml.dom import minidom
import sys

dna = seq1.replace(' ', '')


doc = minidom.Document()

gene = doc.createElement('gene')
gene.setAttribute('acc', acc)
doc.appendChild(gene)

sequence = doc.createElement('sequence')
text = doc.createTextNode(dna)
sequence.appendChild(text)
gene.appendChild(sequence)

for i in range (0, n_cuts):
    enzyme = doc.createElement('enzyme')
    enzyme.setAttribute('name', enz)

    sequence.appendChild(enzyme)

    cut_site = doc.createElement('span')
    text = doc.createTextNode(new_list[i])
    cut_site.appendChild(text)
    enzyme.appendChild(cut_site)

doc.writexml(sys.stdout, addindent='    ', newl='\n')
    
file_handle = open('enz_dummy.xml', 'w')
doc.writexml(file_handle, addindent='   ',newl='\n')
file_handle.close()








    
