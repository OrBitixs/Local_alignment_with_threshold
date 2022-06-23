import numpy as np
from alignment_local_modified import AlignmentLocalModified
from Bio.SubsMat import MatrixInfo

blosum62 = MatrixInfo.blosum62
blosum50 = MatrixInfo.blosum50

seq1 = "INTNINPD"
seq2 = "AINDSQ"

alignment = AlignmentLocalModified(seq1, seq2)
alignment.set_gap(8)
alignment.define_matrix(blosum50)
alignment.set_threshold(3)

alignment.align()
alignment.print_matrix()

print(alignment.last)

alignment.print_traceback_matrix()
print(alignment.last_tuple)

alignment.traceback()
alignment.print_traceback()