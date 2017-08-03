from numpy import mean
from functools import partial
from skbio.alignment import local_pairwise_align_ssw, \
    global_pairwise_align_nucleotide

class DenovoCluster:

    def __init__(self, seqs, similarity_threshold, cluster_fn,
                 aligner=local_pairwise_align_ssw):
        self.seqs = seqs
        self.similarity_threshold = similarity_threshold
        self.cluster_fn = cluster_fn
        self.aligner = aligner

    def build(self):
        clusters = []
        num_alignments = []
        for query_seq in seqs:
            clustered = False
            for i, cluster in enumerate(clusters, start = 1):
                clustered, alignment_results = \
                    self.cluster_fn(query_seq, cluster,
                                    self.similarity_threshold,
                                    self.aligner)
                num_alignments += len(alignment_results)
                if clustered:
                    break
            if clustered:
                for n, s in alignment_results:


def furthest_neighbor(seq, cluster, similarity_threshold, aligner)
    alignment_results = []
    for node in cluster.nodes_iter():
        aln, _, _ = aligner(seq, cluster.node[node]['seq'])
        percent_similarity = 1. - aln[0].distance(aln[1])
        alignment_results.append((node, percent_similarity))
        if percent_similarity < similarity_threshold:
            return False, alignment_results
    return True, alignment_results
