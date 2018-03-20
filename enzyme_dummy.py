#!/usr/bin python3


## dummy file for testing out restriction enzyme cutting sites

import re
from re import finditer

def number_seq(sequence):
    """ Number base pairs in dna seq."""
    count       = 0
    base_count  = {}

    sequence = sequence.replace(' ', '')
    sequence = sequence.replace('\n', '')
    for x in sequence:
        count += 1
        base_count[count] = x
    return base_count
 
        

def GeneSeq(gene):
    """Mapping dna sequence into xml map."""
    if gene in seq:
        seq_map = {'sequence' : seq[gene]}
    return seq_map



def get_enzyme(enzyme, sequence):
    """Get restriction enzyme cutting sequence
        return cutting sites (base numbers) in dna sequence."""
    cut_list = []
    for match in finditer(enzyme, sequence):
        cut_list.append(match.span())
    return cut_list

    

gene_1 = {'AB12345':'ATGGCGGCGCTGTGTCGGACCCGTGCTGTGGCTGCCGAGAGCCATTTTCTGCGAGTGTTTCTCTTCTTCA\
        GGCCCTTTCGGGGTGTAGGCACTGAGAGTGGATCCGAAAGTGGTAGTTCCAATGCCAAGGAGCCTAAGA\
        CGCGCGCAGGCGGTTTCGCGAGCGCGTTGGAGCGGCACTCGGAGCTTCTACAGAAGGTGGAGCCCCTACA\
        GAAGGGTTCTCCAAAAAATGTGGAATCCTTTGCATCTATGCTGAGACATTCTCCTCTTACACAGATGGGA\
        CCTGCAAAGGATAAACTGGTCATTGGACGGATCTTTCATATTGTGGAGAATGATCTGTACATAGATTTTG\
        GTGGAAAGTTTCATTGTGTATGTAGAAGACCAGAAGTGGATGGAGAGAAATACCAGAAAGGAACCAGGGT\
        CCGGTTGCGGCTATTAGATCTTGAACTTACGTCTAGGTTCCTGGGAGCAACAACAGATACAACTGTACTA\
        GAGGCTAATGCAGTTCTCTTGGGAATCCAGGAGAGTAAAGACTCAAGATCGAAAGAAGAACATCATGAAA\
        AAT'}

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

##        

span_list = []
n_cuts = 0
span_list = get_enzyme(fake_cut, seq1)
cut_list = ''
new_list =[]
#print(span_list)
for x in span_list:
    n_cuts += 1
    # have to add one to start and stop because starts numbering at 0
    span   =  str(x[0]+1) + '-' + str(x[1]+1)
    new_list.append(span)
#print(new_list)

#print(n_cuts)



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




## for annotating exon boundaries, etc:


##new_seq = number_seq(seq1)
##
##base = []
##base_num = []
##
##for x in new_seq:
##    base.append(new_seq[x])
##    base_num.append(str(x))
##num_seq = zip(base_num, base)

##
##exon1 = ''
##
##for x in seq1[:20]:
##    exon1 += x
##print(exon1)
##




    
