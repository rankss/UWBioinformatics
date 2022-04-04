import numpy as np
from SequenceAlignment import PWA
from Sequence import Sequence, AASequence, NTSequence
from Score import Score
from Parser import Parser
from Cluster import Cluster
 
def test_SequenceGeneneral():
    NTSeq = Sequence("ACTG", "A")
    assert type(NTSeq) is NTSequence
    assert NTSeq.sequenceType == Sequence.NUCLEOTIDES
    AASeq = Sequence("ACTGQ", "B")
    assert type(AASeq) is AASequence
    assert AASeq.sequenceType == Sequence.AMINO_ACIDS
    AASeq = Sequence("ACGACG", "C", sequenceType=Sequence.AMINO_ACIDS)
    assert type(AASeq) is AASequence
    assert AASeq.sequenceType == Sequence.AMINO_ACIDS
    
def test_SequenceFindSubsequence():
    AASeq = Sequence("ATCCTCGTAATC", "A")
    indices = AASeq.findSubsequence("TC")
    assert indices == [1, 4, 10]

def test_NTSequenceComplement():
    NTSeq = Sequence("ATCGCTAG", "A")
    complementSeq = NTSeq.complement()
    assert complementSeq.sequence == "CTAGCGAT"
    
def test_NTSequenceFindFRSubsequence():
    NTSeq = Sequence("ATCCTCGTAATCGA", "A")
    complementSeq = NTSeq.complement()
    assert complementSeq.sequence == "TCGATTACGAGGAT"
    frIndices = NTSeq.findFRSubsequence("TC")
    assert frIndices["forward"] == [1, 4, 10]
    assert frIndices["reverse"] == [0]
    
def test_NTSequenceTo3AASequences():
    NTSeq = Sequence("ATGAATTAA", "A")
    AASeqs = NTSeq.to3AASequences()
    assert AASeqs[0] == AASequence("MN", "A_translated1")
    assert AASeqs[1] == AASequence("", "A_translated2")
    assert AASeqs[2] == AASequence("EL", "A_translated3")
    
def test_NTSequenceToFR3AASequences():
    pass

# def test_PWAGlobal():
#     matrix = Score(1, -1, 0, -2)
#     alignment = PWA(Sequence("GTCGACGCA", "A"), Sequence("GATTACA", "B"), matrix)
#     alignment.align(PWA.GLOBAL)
#     assert alignment.optimal == -3
#     assert len(alignment.alignments) == 2
#     assert (alignment.dpArray[6] == [-12, -9, -6, -3, -4, -3, 0, -2, -4, -6]).all()

def test_ParserFasta():
    filename = "./Testfiles/test.fasta"
    collection = Parser.Fasta(filename)
    assert collection[0].taxa == "A"
    assert collection[1].taxa == "B"
    assert collection[2].taxa == "C"
    assert collection[0].sequence == "agcaccaaacatgtctaaagctggaaccaaaattactttctttgaagacaaaaactttcaaggccgccactatgacagcgattgcgactgtgcagatttccacatgtacctgagccgctg".upper()
    assert collection[1].sequence == "tgcaccaaacatgtctaaagctggaaccaaaattactttctttgaagacaaaaactttcaaggccgccactatgacagcgattgcgactgtgcagatttccacatgtacctgagccgctg".upper()
    assert collection[2].sequence == "cgcaccaaacatgtctaaagctggaaccaaaattactttctttgaagacaaaaactttcaaggccgccactatgacagcgattgcgactgtgcagatttccacatgtacctgagccgctg".upper()

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
    cluster = Cluster(distances, taxa)
    digraph = cluster.upgma()
    newick = digraph.toNewick()
    assert newick == "(Moth:17.0,(Tuna:14.5,((Turtle:4.0,Chicken:4.0):4.25,(Dog:6.25,(Human:0.5,Monkey:0.5):5.75):2.0):6.25):2.5):0.0"
    
def test_ClusterNJ():
    taxa = ["i", "j", "k", "l"]
    distances = np.array(
        [
            [0, 13, 21, 22],
            [13, 0, 12, 13],
            [21, 12, 0, 13],
            [22, 13, 13, 0]
        ]
    )
    cluster = Cluster(distances, taxa)
    digraph = cluster.nj()
    newick = digraph.toNewick()
    assert newick == "((i:11.0,j:2.0):2.0,(k:6.0,l:7.0):2.0):0.0"
    