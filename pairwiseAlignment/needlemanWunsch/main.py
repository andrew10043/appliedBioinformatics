from NW import *

seqA = 'HEAGAWGHEE'
seqB = 'PAWHEAE'

gapOpen = 8
gapExtend = 8

nw = NW(seqA, seqB, gapOpen, gapExtend)

# Generate F matrix
print(nw.matrix_compute(out="F"))

# Generate T matrix
print(nw.matrix_compute(out="T"))

# Generate final sequence match
print(nw.matrix_compute(out="seq"))

# Generate final matching score
print(nw.matrix_compute(out="score"))

