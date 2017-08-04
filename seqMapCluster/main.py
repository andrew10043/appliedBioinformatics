from skbio import DNA, TabularMSA
import skbio.io
from qiime_default_reference import get_reference_sequences
import random
from skbio.alignment import global_pairwise_align_nucleotide

seqs_16s = []
fraction_to_keep = 0.001
for e in skbio.io.read(get_reference_sequences(), format='fasta',
                       constructor=DNA):
    if random.random() < fraction_to_keep:
        seqs_16s.append(e)

s1 = DNA('AAAAAAAAAA', {'id': 's1'})
s2 = DNA('AAAAATTTTT', {'id': 's2'})
s3 = DNA('AAAAAAACCA', {'id': 's3'})
s4 = DNA('CCCCAATTTT', {'id': 's4'})
s5 = DNA('ACCAAATTTT', {'id': 's5'})
s6 = DNA('AGGAAAAAAA', {'id': 's6'})

aln1 = TabularMSA([s1, s2, s3, s4, s5, s6])

from DenovoCluster import *

cluster = DenovoCluster([s1, s2], 0.7, furthest_neighbor,
                        aligner=global_pairwise_align_nucleotide)

show_clusters(cluster.build())