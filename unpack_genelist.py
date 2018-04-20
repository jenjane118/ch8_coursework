## returns gene identifiers for entire chromosome

# unpacks pickled dictionary object
# prints as xml text file

import pickle, shelve
import Gene_list_class


xmlTemplate = """
<gene>
<geneidentifiers>
    <accession>%(acc)s</accession>
    <genid>%(genid)s</genid>
    <product>%(product)s</product>
    <location>%(location)s</location>
</geneidentifiers>
</gene>"""

f = open('genelist_output.dat', 'rb')
gene_list = pickle.load(f)

print('<?xml version="1.0" encoding="UTF-8"?>\n')
for y in gene_list:
    gene_map = y.Gene_dict()
    print(xmlTemplate%gene_map)


f.close()