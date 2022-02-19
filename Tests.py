from SequenceAlignment import Score, PairwiseAlignment
from Sequence import NUCLEOTIDES
from Matrix import EXAMPLE_MATRIX

def test_pairwiseAlignmentGeneral():
    matrix = Score(1, -1, 0, -1, NUCLEOTIDES)
    print(matrix)
    alignment = PairwiseAlignment("GTCGACGCA", "GATTACA", matrix)
    alignment.Global()
    alignment.Summary()
    
def test_matrixValidation():
    matrix = Score(1, -1, 0, -1, ['A', 'C', 'G'])
    
test_pairwiseAlignmentGeneral()
