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

# Select new query sequences
currentQueries = random.sample(queries, k=10)

# Apply content based heuristics
from CompositionHeuristic import *

comp = CompositionHeuristic(currentQueries, refDB)

# GC content model
print(comp.gc())

# Kmer model (with kmer length of 5)
print(comp.kmer(k=5))