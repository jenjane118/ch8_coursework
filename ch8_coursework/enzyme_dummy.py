#!/usr/bin python3


## dummy file for testing out restriction enzyme cutting sites

import re
import sys



def show_cutDNA(enzyme, sequence):
    """ Get restriction enzyme cutting sequence and return dna seq with indicated
        cutting sites (insert @ symbol)."""
    
    count = 0
    cut = ''
    if enzyme == 'EcoRI':
        cut = 'GAATTC'
    elif enzyme == 'BamHI':
        cut = 'GG.CC'
    elif enzyme == 'BsuMI':
        cut = 'CTCGAG'
    else:
        cut = enzyme
            
    new_seq = sequence.replace(' ', '')
    cut_seq = ''
    p = re.compile(r'('+cut+')')
    it = p.finditer(new_seq)
    for match in it:
        count += 1
    cut_seq = re.sub(r'('+cut+')', r'@\1@', new_seq)
    print(count)
    return cut_seq
    


## You can use either enzyme name ('EcoRI', 'BamHI', or 'BsuMI') or enter custom
## restriction sequence

enz = 'BamHI'

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

# call function with paramenters of enzyme and sequence

new_seq = show_cutDNA(enz, seq1)

# write output in xml to file

from xml.dom import minidom

doc = minidom.Document()

gene = doc.createElement('gene')
gene.setAttribute('acc', acc)
doc.appendChild(gene)

sequence = doc.createElement('sequence')
gene.appendChild(sequence)

enzyme = doc.createElement('enzyme')
enzyme.setAttribute('name', enz)

sequence.appendChild(enzyme)

cut_site = doc.createElement('cut_seq')
text = doc.createTextNode(new_seq)
cut_site.appendChild(text)
enzyme.appendChild(cut_site)

doc.writexml(sys.stdout, addindent='    ', newl='\n')
    
file_handle = open('enz_dummy_out.xml', 'w')
doc.writexml(file_handle, addindent='   ',newl='\n')
file_handle.close()



