import numpy as np
from SequenceAlignment import PairwiseAlignment
from Sequence import Sequence, AASequence, NTSequence
from Score import Score
from Parser import Parser
from Cluster import Cluster, Tree
 
def test_SequenceGeneneral():
    NTSeq = Sequence("ACTG")
    assert type(NTSeq) is NTSequence
    assert NTSeq.sequenceType == Sequence.NUCLEOTIDES
    AASeq = Sequence("ACTGQ")
    assert type(AASeq) is AASequence
    assert AASeq.sequenceType == Sequence.AMINO_ACIDS
    AASeq = Sequence("ACGACG", sequenceType=Sequence.AMINO_ACIDS)
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
    
def test_NTSequenceToAASequence():
    NTSeq = Sequence("ATGAATTAA")
    AASeqs = NTSeq.ToAASequence()
    assert AASeqs[0] == AASequence("MN")
    assert AASeqs[1] == AASequence("")
    assert AASeqs[2] == AASequence("EL")
    print("test_NTSequenceToAASequence: Clear")
    
def test_NTSequenceTo6AASequences():
    pass

def test_pairwiseAlignmentGlobal():
    matrix = Score(1, -1, 0, -2)
    alignment = PairwiseAlignment("GTCGACGCA", "GATTACA", matrix)
    alignment.Align(PairwiseAlignment.GLOBAL)
    assert alignment.optimal == -3
    assert len(alignment.alignments) == 2
    assert (alignment.dpArray[6] == [-12, -9, -6, -3, -4, -3, 0, -2, -4, -6]).all()
    print("test_pairwiseAlignmentGlobal: Clear")

def test_parserFasta():
    filename = "./Testfiles/test.fasta"
    collection = Parser.Fasta(filename)
    assert collection[0].sequenceName == "A"
    assert collection[1].sequenceName == "B"
    assert collection[2].sequenceName == "C"
    assert collection[0] == Sequence("tgcaccaaacatgtctaaagctggaaccaaaattactttctttgaagacaaaaactttcaaggccgccactatgacagcgattgcgactgtgcagatttccacatgtacctgagccgctg", "A")
    assert collection[1] == Sequence("tgcaccaaacatgtctaaagctggaaccaaaattactttctttgaagacaaaaactttcaaggccgccactatgacagcgattgcgactgtgcagatttccacatgtacctgagccgctg", "B")
    assert collection[2] == Sequence("tgcaccaaacatgtctaaagctggaaccaaaattactttctttgaagacaaaaactttcaaggccgccactatgacagcgattgcgactgtgcagatttccacatgtacctgagccgctg", "C")
    print("test_parserFasta: Clear")

def test_ClusterUPGMA():
    # data from table 3 of Fitch and Margoliash, Construction of Phylogenetic trees
    taxa = ["Turtle", "Human", "Tuna", "Chicken", "Moth", "Monkey", "Dog"]
    distances = np.array(
        [
            [0, 19, 27, 8, 33, 18, 13],
            [19, 0, 31, 18, 36, 1, 13],
            [27, 31, 0, 26, 41, 32, 29],
            [8, 18, 26, 0, 31, 17, 14],
            [33, 36, 41, 31, 0, 35, 28],
            [18, 1, 32, 17, 35, 0, 12],
            [13, 13, 29, 14, 28, 12, 0],
        ]
    )
    
    # data from https://en.wikipedia.org/wiki/UPGMA
    taxa = ["a", "b", "c", "d", "e"]
    distances = np.array(
        [
            [0, 17, 21, 31, 23],
            [17, 0, 30, 34, 21],
            [21, 30, 0, 28, 39],
            [31, 34, 28, 0, 43],
            [23, 21, 39, 43, 0]
        ]
    )
    
    cluster = Cluster(distances, taxa)
    newick = Tree.Newick(cluster.UPGMA())
    