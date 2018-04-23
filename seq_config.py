## config file for generating gene list from text file of all gene sequences from chromosome 8

gene_list = []
with open('seq_file.txt', 'r') as f:
    file = f.read().splitlines()
    for x in file:
        gene_list.append(x)


