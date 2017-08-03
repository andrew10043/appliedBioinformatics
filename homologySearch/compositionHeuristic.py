class CompositionHeuristic:

    from skbio.alignment import local_pairwise_align_ssw

    def __init__(self, queries, reference_db, n=5, subset_size=500,
                 aligner=local_pairwise_align_ssw):
        self.queries = queries
        self.reference_db = reference_db
        self.n = n
        self.subset_size = subset_size
        self.aligner = aligner

    # Returns matches based on GC content heuristic
    def gc(self):

        from iab.algorithms import local_alignment_search

        results = []

        # Determine reference database GC content
        reference_db_gc_content = \
            {r.metadata['id']: r.gc_content() for r in self.reference_db}
        for q in self.queries:
            query_gc_content = q.gc_content()
            database_subset = []
            for r in self.reference_db:
                ref_gc_content = reference_db_gc_content[r.metadata['id']]
                database_subset.append((abs(ref_gc_content - query_gc_content),
                                        r))
            database_subset.sort(key=lambda x: x[0])
            database_subset = \
                [e[1] for e in database_subset[:self.subset_size]]
            results.append(local_alignment_search(
                [q], database_subset, n=self.n, aligner=self.aligner))

        import pandas as pd

        return print(pd.concat(results))

    # Returns matches based on kmer heuristic
    def kmer(self, k=5):

        from iab.algorithms import local_alignment_search
        from iab.algorithms import fraction_shared_kmers

        results = []

        # Determine reference kmer frequency
        reference_db_kmer_frequencies = \
            {r.metadata['id']:
                 r.kmer_frequencies(k=k, overlap=True)
             for r in self.reference_db}
        for q in self.queries:
            query_kmer_frequency = q.kmer_frequencies(k=k, overlap=True)
            database_subset = []
            for r in self.reference_db:
                ref_kmer_frequency = \
                    reference_db_kmer_frequencies[r.metadata['id']]
                s = fraction_shared_kmers(query_kmer_frequency,
                                          ref_kmer_frequency)
                database_subset.append((s, r))
            database_subset.sort(key=lambda x: x[0], reverse = True)
            database_subset = \
                [e[1] for e in database_subset[:self.subset_size]]
            results.append(local_alignment_search(
                [q], database_subset, n=self.n, aligner=self.aligner))

        import pandas as pd

        return print(pd.concat(results))



