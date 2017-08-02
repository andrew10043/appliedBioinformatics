# Import BLOSUM50 amino acid substitution matrix
from iab.algorithms import blosum50

import numpy as np

import pandas as pd

# Generate amino acid sequences
seqA = "HEAGAWGHEE"
seqB = "PAWHEAE"

# Create the F matrix with zeros at all positions
numRows = len(seqB) + 1
numCols = len(seqA) + 1

rowLabs = list(seqB)
rowLabs.insert(0, "")

colLabs = list(seqA)
colLabs.insert(0, "")

emptyF = np.zeros(shape=(numRows, numCols), dtype=np.int)
F = pd.DataFrame(emptyF, columns=colLabs, index=rowLabs)

# Initialize F using the following rules:
# F(0, 0) = 0
# F(i, 0) = F(i - 1, 0) - d
# F(0, j) = F(0, j - 1) - d
# Where d is a penalty for insertion of a gap character to align the seqs
d = 8

for i in range(1, F.shape[0]):
    F.iloc[i, 0] = F.iloc[i - 1, 0] - d

for j in range(1, F.shape[1]):
    F.iloc[0, j] = F.iloc[0, j - 1] - d

# Define characters for T table
trace = {'diag': '\u2196',
         'left': '\u2190',
         'up': '\u2191'
         }

startDot = '\u25cf'

# Generate T table
emptyT = np.zeros(shape=(numRows, numCols), dtype=np.int)
T = pd.DataFrame(emptyT, columns=colLabs, index=rowLabs)

# Intialize T table such that each cell contains an arrow pointing towards
# the cell whose F score it depends on
T.iloc[0, 0] = startDot
for i in range(1, T.shape[0]):
    T.iloc[i, 0] = upArrow
for j in range(1, T.shape[1]):
    T.iloc[0, j] = leftArrow

# Fill the F matrix with the following rules:
# F(i, j) = max of:
# 1. F(i - 1, j - 1) + s(ci, cj) where s(ci, cj) is the substitution score
# 2. F(i - 1, j) - d
# 3. F(i, j - 1) - d

# Concurrently fill the T matrix based on which score defines the current cell
# of the F matrix

for i in range(1, F.shape[0]):
    for j in range(1, F.shape[1]):
        subPenalty = blosum50[F.index[i]][F.columns[j]]
        diagScore = F.iloc[i - 1, j - 1] + subPenalty
        leftScore = F.iloc[i, j - 1] - d
        upScore = F.iloc[i - 1, j] - d
        maxScore = max(diagScore, leftScore, upScore)
        F.iloc[i, j] = maxScore
        if maxScore == diagScore:
            T.iloc[i, j] = trace['diag']
        elif maxScore == upScore:
            T.iloc[i, j] = trace['up']
        elif maxScore == leftScore:
            T.iloc[i, j] = trace['left']

