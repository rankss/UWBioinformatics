from Constants import *
from Logger import Logger
from Sequence import *
from Score import Score

class PairwiseAlignment:
    """[summary]
    """
    def __init__(self, sequence1: Sequence, sequence2: Sequence, score: Score):
        self.logger = Logger()
        self.hsequence = sequence1.sequence
        self.vsequence = sequence2.sequence
        self.score = score
        self.alignment = None
        return

    def __str__(self):
        # Prints alignment
        return ""

    def Global(self):
        """Performs Needleman-Wunsch for sequence pair
        """
        hlen = len(self.hsequence)
        vlen = len(self.vsequence)

        return

    def Local(self):
        """Performs Smith-Waterman for sequence pair
        """

        return

score = Score(MatrixConstants.EXAMPLE_MATRIX, BlockConstants.NUCLEOTIDES, 2, 1)
DNASeq1 = DNASequence("CAAGAC")
DNASeq2 = DNASequence("GAAC")
alignment = PairwiseAlignment(DNASeq1, DNASeq2, score)
alignment.Global()