from SequenceAlignment import Score, PairwiseAlignment
from Constants import *

def test_pairwiseAlignmentGeneral():
    matrix = Score(Matrix.EXAMPLE_MATRIX, Block.NUCLEOTIDES, 0, -1)
    alignment = PairwiseAlignment("GTCGACGCA", "GATTACA", matrix)
    alignment.Local()
    alignment.Summary()

test_pairwiseAlignmentGeneral()
