import seaborn as sns
from functools import partial
import numpy as np
import random
import skbio
from skbio import DistanceMatrix
from skbio.sequence.distance import hamming
from skbio.alignment import global_pairwise_align_nucleotide
import ete3
from iab.algorithms import evolve_generations, random_sequence, \
    progressive_msa, kmer_distance


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

ts = ete3.TreeStyle()
ts.show_leaf_name = True
ts.scale = 250
ts.branch_vertical_margin = 15

t = ete3.Tree()
t.populate(10)

# Plot distance matrices using various metrics for distance

# Alignment-Free (i.e. kmer)
kmer_dm = DistanceMatrix.from_iterable(sequences, metric=kmer_distance,
                                       key='id')
plot = kmer_dm.plot(cmap='Greens', title='3mer distances between sequences')
plot.savefig('kmer_distance_matrix.png')

# Aligment-based (Hamming distance) - need to align seqs first

hamming_dm = DistanceMatrix.from_iterable(sequences_aligned, metric=hamming,
                                          key='id')
plot = hamming_dm.plot(cmap='Greens',
                       title='Hamming distances between sequences')
plot.savefig('hamming_distance_matrix.png')


# Jukes-Cantor correction of observed distances between sequences


def jc_correction(p):
    # Return max of original value and JC correction
    # This avoids the issues of returning 'nan' values if 4/3p > 1
    return max(p, (-3/4) * np.log(1 - (4*p/3)))


def jc_correct_dm(dm):
    result = np.zeros(dm.shape)
    for i in range(dm.shape[0]):
        for j in range(i):
            result[i, j] = result[j, i] = jc_correction(dm[i, j])
    return skbio.DistanceMatrix(result, ids=dm.ids)


jc_dm = jc_correct_dm(hamming_dm)

plot = jc_dm.plot(cmap='Greens',
                  title='JC-corrected Hamming Distances between sequences')
plot.savefig('jc_corrected_distance_matrix.png')

