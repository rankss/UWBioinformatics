from SequenceAlignment import Score, PairwiseAlignment
from Sequence import NUCLEOTIDES
from Model import Model
import numpy as np

def test_pairwiseAlignmentGeneral():
    matrix = Score(1, -1, 0, -1, NUCLEOTIDES)
    print(matrix)
    alignment = PairwiseAlignment("GTCGACGCA", "GATTACA", matrix)
    alignment.Global()
    alignment.Summary()
    
def test_matrixValidation():
    try:
        matrix = Score(1, -1, 0, -1, ['A', 'C', 'G'])
    except:
        print(f"Error thrown")
    
def test_phasePortrait():
    model = Model()
    
    def dX(x, u):
        return u*x - x**3
    
    def dY(y):
        return -y
    
    X, Y = np.meshgrid(np.arange(-3, 3, 0.01), np.arange(-3, 3, 0.01))
    model.PhasePortrait(dX, dY, (X, -1), (Y,), X, Y)
    return
        
test_pairwiseAlignmentGeneral()
