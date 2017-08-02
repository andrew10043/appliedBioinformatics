from SW import *

seqA = "HEAGAWGHEE"
seqB = "PAWHEAE"

gapOpen = 8
gapExtend = 8

sw = SW(seqA, seqB, gapOpen, gapExtend)

# Generate the F matrix
print(sw.matrix_compute(out="F"))

# Generate the T matrix
print(sw.matrix_compute(out="T"))

# Generate the final sequence match
print(sw.matrix_compute(out="seq"))

# Generate the matching score
print(sw.matrix_compute(out="score"))