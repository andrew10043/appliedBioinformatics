from iab.algorithms import blosum50
import numpy as np
import pandas as pd


class NW:

    def __init__(self, seqA, seqB, gapOpen, gapExtend):
        self.seqA = seqA
        self.seqB = seqB
        self.gapOpen = gapOpen
        self.gapExtend = gapExtend

    # Returns F and T matrix, sequence alignment and final score
    def matrix_compute(self, out):
        numRows = len(self.seqB) + 1
        numCols = len(self.seqA) + 1

        rowLabs = list(self.seqB)
        rowLabs.insert(0, "")

        colLabs = list(self.seqA)
        colLabs.insert(0, "")

        emptyF = np.zeros(shape=(numRows, numCols), dtype=np.int)
        F = pd.DataFrame(emptyF, columns=colLabs, index=rowLabs)

        # Initialize F using the following rules:
        # F(0, 0) = 0
        # F(i, 0) = F(i - 1, 0) - d
        # F(0, j) = F(0, j - 1) - d
        # Where d is gapOpen penalty
        d = self.gapOpen

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

        # Initialize T  such that each cell contains an arrow pointing towards
        # the cell whose F score it depends on
        T.iloc[0, 0] = startDot
        for i in range(1, T.shape[0]):
            T.iloc[i, 0] = trace['up']
        for j in range(1, T.shape[1]):
            T.iloc[0, j] = trace['left']

        # Fill the F matrix with the following rules:
        # F(i, j) = max of:
        # 1. F(i - 1, j - 1) + s(ci, cj); s(ci, cj) is the substitution score
        # 2. F(i - 1, j) - d
        # 3. F(i, j - 1) - d
        # Gap extend (e) used in place of d if more than one adjacent gap

        # Fill the T matrix based on which score defines the current cell
        # of the F matrix

        # Default in the case of ties is to take the left score, then
        # diag, then up (this is the default in the source code) - changes to
        # this order will change the final sequence match.

        d = self.gapOpen
        e = self.gapExtend

        for i in range(1, F.shape[0]):
            for j in range(1, F.shape[1]):

                subPenalty = blosum50[F.index[i]][F.columns[j]]
                diagScore = F.iloc[i - 1, j - 1] + subPenalty

                if T.iloc[i, j - 1] == trace['up']:
                    leftScore = F.iloc[i, j - 1] - e
                else:
                    leftScore = F.iloc[i, j - 1] - d

                if T.iloc[i - 1, j] == trace['left']:
                    upScore = F.iloc[i - 1, j] - e
                else:
                    upScore = F.iloc[i - 1, j] - d

                maxScore = max(diagScore, leftScore, upScore)

                F.iloc[i, j] = maxScore

                if maxScore == leftScore:
                    T.iloc[i, j] = trace['left']
                elif maxScore == diagScore:
                    T.iloc[i, j] = trace['diag']
                elif maxScore == upScore:
                    T.iloc[i, j] = trace['up']

        i = T.shape[0] - 1
        j = T.shape[1] - 1

        alignB = []
        alignA = []

        while (i, j) != (0, 0):
            if T.iloc[i, j] == trace['diag']:
                alignB.insert(0, T.index[i])
                alignA.insert(0, T.columns[j])
                i = i - 1
                j = j - 1
            elif T.iloc[i, j] == trace['up']:
                alignB.insert(0, T.index[i])
                alignA.insert(0, "-")
                i = i - 1
                j = j
            elif T.iloc[i, j] == trace['left']:
                alignB.insert(0, "-")
                alignA.insert(0, T.columns[j])
                i = i
                j = j - 1

        finalSeq = np.matrix([["".join(alignA)], ["".join(alignB)]])

        finalScore = F.iloc[F.shape[0] - 1, F.shape[1] - 1]

        if out == "F":
            return F
        elif out == "T":
            return T
        elif out == "seq":
            return finalSeq
        elif out == "score":
            return finalScore