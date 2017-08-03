import numpy as np
import seaborn as sns
import random
from iab.algorithms import evolve_generations
from iab.algorithms import random_sequence
from iab.algorithms import progressive_msa
from skbio.alignment import global_pairwise_align_nucleotide
import skbio

# Generate random last common ancestor sequence
sequence = random_sequence(skbio.DNA, 50)

# Simulate 10 generations
# We set indel probability to 0 here so we won't have to align the seqs
# before constructing a tree (to save run time for the example)
indel_probability = 0.0

sequences = evolve_generations(sequence, generations=10,
                               substitution_probability=0.1,
                               indel_probability=0.0,
                               increased_rate_probability=0.1,
                               verbose=False)

sequences = random.sample(sequences, 25)

# Screening step to align sequences if indel prob > 0
if indel_probability == 0:
    sequences_aligned = sequences
else:
    gpa = partial(global_pairwise_align_nucleotide,
                  penalize_terminal_gaps=True)
    sequences_aligned = progressive_msa(sequences,
                                        pairwise_aligner=gpa)

import ete3
ts = ete3.TreeStyle()
ts.show_leaf_name = True
ts.scale = 250
ts.branch_vertical_margin = 15

t = ete3.Tree()
t.populate(10)
t.render("test.png")

# Plot distance matrices using various metrics for distance

# Alignment-Free (i.e. kmer)
from iab.algorithms import kmer_distance
from skbio import DistanceMatrix

kmer_dm = DistanceMatrix.from_iterable(sequences, metric=kmer_distance,
                                       key='id')
plot = kmer_dm.plot(cmap='Greens', title='3mer distances between sequences')
plot.savefig('kmer_distance_matrix.png')

# Aligment-based (Hamming distance) - need to align seqs first
from skbio.sequence.distance import hamming
hamming_dm = DistanceMatrix.from_iterable(sequences_aligned, metric=hamming,
                                          key='id')
plot = hamming_dm.plot(cmap='Greens',
                       title='Hamming distances between sequences')
plot.savefig('hamming_distance_matrix.png')


