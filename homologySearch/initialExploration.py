from skbio.alignment import local_pairwise_align_ssw

import qiime_default_reference as qdr
import skbio

# Load taxonomic data
refTax = {}
for e in open(qdr.get_reference_taxonomy()):
    seqID, seqTax = e.strip().split('\t')
    refTax[seqID] = seqTax

# Load reference sequences and associate with taxonomic data as metadata
refDB = []
for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta',
                       constructor=skbio.DNA):
    seqTax = refTax[e.metadata['id']]
    e.metadata['taxonomy'] = seqTax
    refDB.append(e)

# Utilize random selection of reference sequences for examples
import random
refDB = random.sample(refDB, k=5000)

# Select a random subset of the reference sequences to use as mock queries
# No annotation with taxonomic metadata will occur
# Sequences will also be trimmed to simulate partial gene sequences
queries =[]
for e in skbio.io.read(qdr.get_reference_sequences(), format='fasta',
                       constructor=skbio.DNA):
    e = e[100:300]
    queries.append(e)
queries = random.sample(queries, k=500)

# Import local_alignment_search
# This function utilizes local pairwise alignment (i.e. the SW algorithm)
from iab.algorithms import local_alignment_search

# Run local_alignment_search for 4 randomly selected queries
import time
startTime = time.time()
currentQueries = random.sample(queries, k=4)
results = local_alignment_search(currentQueries, refDB)
stopTime = time.time()
print("Run time: " + str((stopTime - startTime)/len(currentQueries)) +
      " seconds per query.")
print(results)

# Print actual taxonomic metadata for the query sequences to compare
for q in currentQueries:
    qID = q.metadata['id']
    print("Known taxonomy for query " + str(qID) + ":\n" +
          str(refTax[qID]))
