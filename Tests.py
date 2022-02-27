from SequenceAlignment import Score, PairwiseAlignment
from Sequence import Sequence, AASequence, NTSequence
import numpy as np

def test_pairwiseAlignmentGLOBAL():
    matrix = Score(1, -1, 0, -2)
    alignment = PairwiseAlignment("GTCGACGCA", "GATTACA", matrix)
    alignment.Align(PairwiseAlignment.GLOBAL)
    assert alignment.optimal == -3
    assert len(alignment.alignments) == 2
    assert (alignment.dpArray[6] == np.array([-12, -9, -6, -3, -4, -3, 0, -2, -4, -6])).all()
    print("test_pairwiseAlignmentGLOBAL: Clear")
    
def test_SequenceGeneneral():
    ntseq = Sequence("ACTG")
    assert type(ntseq) is NTSequence
    assert ntseq.sequenceType == Sequence.NUCLEOTIDES
    aaseq = Sequence("ACTGQ")
    assert type(aaseq) is AASequence
    assert aaseq.sequenceType == Sequence.AMINO_ACIDS
    print("test_SequenceGeneneral: Clear")
