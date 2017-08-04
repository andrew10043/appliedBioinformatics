from skbio.alignment import local_pairwise_align_ssw
import networkx as nx

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
                    cluster.add_node(query_seq.metadata['id'], seq=query_seq)
                    cluster.add_edge(query_seq.metadata['id'], n,
                                     percent_similarity=s)
                    cluster.graph['node-order'].append(query_seq.metadata['id'])
            else:
                new_cluster = nx.Graph(id="OTU %d" % (len(clusters) + 1))
                new_cluster.add_node(query_seq.metadata['id'],
                                     seq=query_seq)
                new_cluster.graph['node-order'] = [query_seq.metadata['id']]
                clusters.append(new_cluster)
        return clusters, num_alignments


def furthest_neighbor(seq, cluster, similarity_threshold, aligner):
    alignment_results = []
    for node in cluster.nodes_iter():
        aln, _, _ = aligner(seq, cluster.node[node]['seq'])
        percent_similarity = 1. - aln[0].distance(aln[1])
        alignment_results.append((node, percent_similarity))
        if percent_similarity < similarity_threshold:
            return False, alignment_results
    return True, alignment_results

def show_clusters(clusters, print_clusters=True, plot_clusters=True,
                  plot_labels=False):
    G = nx.Graph()
    for c in clusters:
        G = nx.union(G, c)
        if print_clusters: print("%s: %s" %
                                 (c.graph['id'],
                                  [s for s in c.graph['node-order']]))
    if plot_clusters:
        pos=nx.spring_layout(G)
        nx.draw_networkx_nodes(G, pos, node_color='w')
        if G.number_of_edges() > 0:
            nx.draw_networkx_edges(G, pos)
        if plot_labels:
            nx.draw_networkx_labels(G, pos,
                                    labels={d:G.node[d]['seq'].metadata['id']
                                            for d in G.nodes()})
        _ = plt.axis('off')