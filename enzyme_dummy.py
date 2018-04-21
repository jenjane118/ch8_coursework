#!/usr/bin python3


## dummy file for testing out restriction enzyme cutting sites

import re
import sys
import seq_module



from xml.dom import minidom

def show_cutDNA(acc, sequence, enzyme=None):
    """ Get restriction enzyme cutting sequence and return dna seq with indicated
        cutting sites (insert @ symbol)."""

    enz_dict    = {
        'EcoRI': 'GAATTC', 'BamHI': 'GGATCC',
        'BsuMI': 'CTCGAG', 'HindIII': 'AAGCTT',
        'EcoRV': 'GATATC'}

    cut_dict    = {}

    ## if a custom cleavage site is included, will indicate number and position of cleavage sites in sequence
    if enzyme   != None:
        count   = 0
        cut     = enzyme
        new_seq = sequence.replace(' ', '')
        p       = re.compile(r'(' + cut + ')')
        it      = p.finditer(new_seq)
        for match in it:
            count += 1
        cut_seq             = re.sub(r'(' + cut + ')', r'@\1@', new_seq)
        cut_dict[enzyme]    = (count, cut_seq)

    ## for each enzyme in dictionary, it will search and return only the enzymes which have cleavage sites in sequence
    for enzyme in enz_dict:
        count = 0
        cut = enz_dict[enzyme]
        new_seq = sequence.replace(' ', '')
        p = re.compile(r'('+cut+')')
        it = p.finditer(new_seq)
        for match in it:
            count += 1
        cut_seq = re.sub(r'('+cut+')', r'@\1@', new_seq)
        if count != 0:
            cut_dict[enzyme] = (count, cut_seq)
    return cut_dict

## Main



## You can use either enzyme name ('EcoRI', 'BamHI', or 'BsuMI') or enter custom
## restriction sequence


#enz = 'BamHI'

custom = 'GTGGCG'


acc  = 'AB12345'
seq1 = "ATGGCAATGCAGTGGCGCTGTGTCGGACCCGTGCTGTGGCTGCCGAGAGCCATTTTCTGCGAGTGTTTCTCTTCTTCAGGCCCTTTCGGGGTGTAGGCACTGAGAGTGGATCCGAAAGTGGTAGTTCCAATGCCAAGGAGCCTAAGACGCGCGCAGGCGGTTTCGCGAGCGCGTTGGAGCGGCACTCGGAGCTTCTACAGAAGGTGGAGCCCCTACAGAAGGGTTCTCCAAAAAATGTGGAATCCTTTGCATCTATGCTGAGACATTCTCCTCTTACACAGATGGGACCTGCAAAGGATAAACTGGTCATTGGACGGATCTTTCATATTGTGGAGAATGATCTGTACATAGATTTTGGTGGAAAGTTTCATTGTGTATGTAGAAGACCAGAAGTGGATGGAGAGAAATACCAGAAAGGAACCAGGGTCCGGTTGCGGCTATTAGATCTTGAACTTACGTCTAGGTTCCTGGGAGCAACAACAGATACAACTGTACTAGAGGCTAATGCAGTTCTCTTGGGAATCCAGGAGAGTAAAGACTCAAGATCGAAAGAAGAACATCATGAAAAAT"

exons = [('AB12345', 55, 257)]


## call function with paramenters of enzyme and sequence, show whether enzymes cut inside coding region
exon_seq    = seq_module.annotateSeq(acc, seq1, exons)
exon_cut    = show_cutDNA(acc, exon_seq, custom)

#for k,v in exon_cut.items():
#    print(k, ':', v[0], '\n', v[1])




## write output in xml to file

doc = minidom.Document()

seq_data = doc.createElement('seq_data')
doc.appendChild(seq_data)

for k,v in exon_cut.items():
    enz = k
    no = str(v[0])
    new_seq = v[1]

    gene = doc.createElement('gene')
    gene.setAttribute('acc', acc)

    seq_data.appendChild(gene)

    sequence = doc.createElement('sequence')
    gene.appendChild(sequence)

    enzyme = doc.createElement('enzyme')
    enzyme.setAttribute('name', enz)
    enzyme.setAttribute('number', no)

    sequence.appendChild(enzyme)

    cut_site = doc.createElement('cut_seq')
    text = doc.createTextNode(new_seq)
    cut_site.appendChild(text)
    enzyme.appendChild(cut_site)

    #doc.writexml(sys.stdout, addindent='    ', newl='\n')

file_handle = open('enz_dummy_out.xml', 'w')
doc.writexml(file_handle, addindent='   ',newl='\n')
file_handle.close()



