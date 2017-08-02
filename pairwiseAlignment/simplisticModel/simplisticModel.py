# Example of an overly simplistic pairwise sequence alignment algorithm

import numpy as np

import pandas as pd

# Define sequences
seqA = "ACCGGTGACCTAAC"
seqB = "ACCGGTCCGTAAC"

# Generate comparative data frame; initialize with all zeros
numRows = len(seqB)
numCols = len(seqA)
empty = np.zeros(shape=(numRows, numCols), dtype=np.int)
data = pd.DataFrame(empty, columns=list(seqA), index=list(seqB))

# Replace matching entries with a "1"
for colLab, col in data.iteritems():
    for rowLab, row in data.iterrows():
        if colLab == rowLab:
            data.loc[rowLab, colLab] = 1
        elif colLab != rowLab:
            data.loc[rowLab, colLab] = 0

summedData = data.copy()

# For each cell (starting in row 2, col 2), add the value of the cell to the
# left upward diagonal. This allows us to document local areas of sequence
# matches.
for i in range(1, summedData.shape[0]):
    for j in range(1, summedData.shape[1]):
        if summedData.iloc[i, j] > 0:
            summedData.iloc[i, j] += summedData.iloc[i-1, j-1]

# Print final data frame, with zeros replaced by blanks for ease of viewing
print(summedData.replace(0, ""))

# Identify the longest diagonal (matching sequence)
print("The longest diagonal is " + str(summedData.values.max()) +
      " characters long.")

