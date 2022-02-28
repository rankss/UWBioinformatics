from SequenceAlignment import Score, PairwiseAlignment
from Sequence import Sequence, AASequence, NTSequence
import numpy as np
 
def test_SequenceGeneneral():
    NTSeq = Sequence("ACTG")
    assert type(NTSeq) is NTSequence
    assert NTSeq.sequenceType == Sequence.NUCLEOTIDES
    AASeq = Sequence("ACTGQ")
    assert type(AASeq) is AASequence
    assert AASeq.sequenceType == Sequence.AMINO_ACIDS
    AASeq = Sequence("ACGACG", Sequence.AMINO_ACIDS)
    assert type(AASeq) is AASequence
    assert AASeq.sequenceType == Sequence.AMINO_ACIDS
    print("test_SequenceGeneneral: Clear")
    
def test_SequenceFindSubsequence():
    AASeq = Sequence("ATCCTCGTAATC")
    indices = AASeq.FindSubsequence("TC")
    assert indices == [1, 4, 10]
    print("test_SequenceFindSubsequence: Clear")

def test_NTSequenceComplement():
    NTSeq = Sequence("ATCGCTAG")
    complementSeq = NTSeq.Complement()
    assert complementSeq.sequence == "CTAGCGAT"
    print("test_NTSequenceComplement: Clear")
    
def test_NTSequenceFindFRSubsequence():
    NTSeq = Sequence("ATCCTCGTAATCGA")
    complementSeq = NTSeq.Complement()
    assert complementSeq.sequence == "TCGATTACGAGGAT"
    frIndices = NTSeq.FindFRSubsequence("TC")
    assert frIndices["forward"] == [1, 4, 10]
    assert frIndices["reverse"] == [0]
    print("test_NTSequenceFindFRSubsequence: Clear")

def test_pairwiseAlignmentGlobal():
    matrix = Score(1, -1, 0, -2)
    alignment = PairwiseAlignment("GTCGACGCA", "GATTACA", matrix)
    alignment.Align(PairwiseAlignment.GLOBAL)
    assert alignment.optimal == -3
    assert len(alignment.alignments) == 2
    assert (alignment.dpArray[6] == np.array([-12, -9, -6, -3, -4, -3, 0, -2, -4, -6])).all()
    print("test_pairwiseAlignmentGlobal: Clear")
