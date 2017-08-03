import collections

from skbio.alignment import local_pairwise_align_ssw

import pandas as pd
import numpy as np

def evaluate_search(queries, refDB, refTax, searchFunction, taxLevel,
                    n=5, aligner=local_pairwise_align_ssw):
    startTime = time.time()
    searchResults = searchFunction(queries, refDB, n=n, aligner=aligner)
    stopTime = time.time()
    runTime = stopTime - startTime
    perQueryRuntime = runTime/len(queries)
    data = []
    indices = []
    for q in queries:
        qID = q.metadata['id']
        indices.append(qID)
        qKnownTax = tuple(refTax[qID].split('; ')[:taxLevel])
        qObsTaxAll = collections.Counter()
        for e in searchResults['reference taxonomy'][qID]:
            qObsTaxAll[tuple(e.split('; ')[:taxLevel])] += 1
        qObsTaxOne = qObsTaxAll.most_common()[0][0]
        data.append((qKnownTax, qObsTaxOne))
    index = pd.Index(indices, name='Query ID')
    data = pd.DataFrame(data, index=index,
                        columns=['Known taxonomy', 'Observed taxonomy'])
    numCorrect = np.sum(data['Known taxonomy'] == data['Observed taxonomy'])
    fracCorrect = numCorrect / data.shape[0]
    return perQueryRuntime, fracCorrect, data



