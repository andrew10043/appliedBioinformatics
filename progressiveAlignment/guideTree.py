from skbio import DNA
from iab.algorithms import kmer_distance

# Build initial guide distance matrix based on kmer metric
query_sequences = [DNA("ACCGGTGACCAGTTGACCAGT", {"id": "s1"}),
                   DNA("ATCGGTACCGGTAGAAGT", {"id": "s2"}),
                   DNA("GGTACCAAATAGAA", {"id": "s3"}),
                   DNA("GGCACCAAACAGAA", {"id": "s4"}),
                   DNA("GGCCCACTGAT", {"id": "s5"})]

from skbio import DistanceMatrix

guide_dm = DistanceMatrix.from_iterable(query_sequences,
                                        metric=kmer_distance, key='id')

# Cluster the sequences using the unweighted pair group method with arithmetic
# mean (UPGMA) to create a rooted dendrogram

from scipy.cluster.hierarchy import average, dendrogram, to_tree

for q in query_sequences:
    print(q)

guide_lm = average(guide_dm.condensed_form())
guide_d = dendrogram(guide_lm, labels=guide_dm.ids, orientation = 'right',
                     link_color_func=lambda x: 'black')
guide_tree = to_tree(guide_lm)

# Alternatively, we can use a pre-defined function

from iab.algorithms import guide_tree_from_sequences

t = guide_tree_from_sequences(query_sequences, display_tree=False)