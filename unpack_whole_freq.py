## returns codon usage info for entire chromosome

# unpickling dictionaries generated from whole_genome_freq.py

import pickle, shelve

f = open('whole_genome_usage.dat', 'rb')
usage_ratio = pickle.load(f)
usage_percent = pickle.load(f)

print (usage_ratio)
print(usage_percent)

f.close()

