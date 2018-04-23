## config file for generating gene list from text file of all gene sequences from chromosome 8

gene_list = []
with open('seq_file.txt', 'r') as f:
    file = f.read().splitlines()
    for x in file:
        gene_list.append(x)


# with cnx.cursor() as cursor:
#     query4 = " SELECT accession,  positions, translation FROM coding_sequence;"
#     cursor.execute(query4)
#     coding_seq = cursor.fetchall()
#     print(coding_seq)  # sanity check
#     data = cursor.fetchall()
#     for row in data:
        #do
        #stuff
## make into dictionary? is this necessary?