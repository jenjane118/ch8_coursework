import seq_module


# for each gene in gene_list:
    # run codingSeq to get coding sequence
    # run codon freq (output freq table/dict)
    # update values in whole genome dictionary



## for calculating whole genome codon frequency
genome_freq = {}
for gene in gene_list:
    seq         = getSequence(gene)
    code_seq    = codingSeq(seq)
    code_freq   = codonFreq
    for key in code_freq:
        if key in genome_freq:
            genome_freq[key] += code_freq[key]
        else:
            genome_freq = {key : code_freq[key]}
for k, v in genome_freq.items():
    print(k, v)

