from iab.algorithms import blosum50
import numpy as np
import pandas as pd

class SW:

    def __init__(self, seqA, seqB, gapOpen, gapExtend):
        self.seqA = seqA
        self.seqB = seqB
        self.gapOpen = gapOpen
        self.gapExtend = gapExtend

    def matrix_compute(self, out):
        numRows = len(self.seqB) + 1
        numCols = len(self.seqA) + 1

        rowLabs = list(self.seqB)
        rowLabs.insert(0, "")

        colLabs = list(self.seqA)
        colLabs.insert(0, "")

        # In the SW local matching algorithm, F is initialized with
        # all zeros (unlike NW)
        emptyF = np.zeros(shape=(numRows, numCols), dtype=np.int)
        F = pd.DataFrame(emptyF, columns=colLabs, index=rowLabs)

        # Define characters for T table
        trace = {'diag': '\u2196',
                 'left': '\u2190',
                 'up': '\u2191'
                 }

        startDot = '\u25cf'

        # Build and initialize the T matrix
        # Given the F matrix initialization with all zeros, the T matrix
        # is inialized with startDots in the entire first row and column

        emptyT = np.zeros(shape=(numRows, numCols), dtype=np.int)
        T = pd.DataFrame(emptyT, columns=colLabs, index=rowLabs)

        for i in range(0, T.shape[0]):
            T.iloc[i, 0] = startDot
        for j in range(0, T.shape[1]):
            T.iloc[0, j] = startDot

        # Fill the F matrix with the following rules:
        # F(i, j) = max of:
        # 1. 0
        # 2. F(i - 1, j - 1) + s(ci, cj); s(ci, cj) is the substitution score
        # 3. F(i - 1, j) - d
        # 4. F(i, j - 1) - d
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

                maxScore = max(diagScore, leftScore, upScore, 0)

                F.iloc[i, j] = maxScore

                if maxScore == 0:
                    T.iloc[i, j] = startDot
                elif maxScore == leftScore:
                    T.iloc[i, j] = trace['left']
                elif maxScore == diagScore:
                    T.iloc[i, j] = trace['diag']
                elif maxScore == upScore:
                    T.iloc[i, j] = trace['up']

        # Traceback in SW starts with the maximum F cell value, unlike
        # NW which always starts in the bottom right corner
        maxValue = 0.0
        max_i = 0
        max_j = 0
        for i in range(F.shape[0]):
            for j in range(F.shape[1]):
                if F.iloc[i, j] > maxValue:
                    max_i, max_j = i, j
                    maxValue = F.iloc[i, j]

        # In SW, traceback terminates when it reaches a startDot, which
        # unlike NW can occur in locations other than the upper right corner

        alignB = []
        alignA = []

        i = max_i
        j = max_j

        while T.iloc[i, j] != startDot:
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

        finalScore = F.iloc[max_i, max_j]

        if out == "F":
            return F
        elif out == "T":
            return T
        elif out == "seq":
            return finalSeq
        elif out == "score":
            return finalScore
