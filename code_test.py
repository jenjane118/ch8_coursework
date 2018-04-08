with open ('seq_file1.txt', 'r') as f:
    file = f.read().splitlines()


CodonsDict = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
        'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}


codon_list = []
codon = ''
#print(file)
for x in file:

    #print(x)
    codon += x
    print(codon)
    if len(codon) == 3:
        print(codon)
        # codon_list.append(x[i]:x[i+3])
        codon_list.append(codon)
        codon = ''
    else:
        continue
#print(codon_list)
# for x in codon_list:
#     #count = 0
#     if x in CodonsDict:
#         CodonsDict[x] += 1
# return CodonsDict

SynCodons = {
        'C' : ['TGT', 'TGC'],
        'D' : ['GAT', 'GAC'],
        'S' : ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'Q' : ['CAA', 'CAG'],
        'M' : ['ATG'],
        'N' : ['AAC', 'AAT'],
        'P' : ['CCT', 'CCG', 'CCA', 'CCC'],
        'K' : ['AAG', 'AAA'],
        'T' : ['ACC', 'ACA', 'ACG', 'ACT'],
        'F' : ['TTT', 'TTC'],
        'A' : ['GCA', 'GCC', 'GCG', 'GCT'],
        'G' : ['GGT', 'GGG', 'GGA', 'GGC'],
        'I' : ['ATC', 'ATA', 'ATT'],
        'L' : ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
        'H' : ['CAT', 'CAC'],
        'R' : ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'W' : ['TGG'],
        'V' : ['GTA', 'GTC', 'GTG', 'GTT'],
        'E' : ['GAG', 'GAA'],
        'Y' : ['TAT', 'TAC'],
        'STOP': ['TAG', 'TGA', 'TAA']}

codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}
