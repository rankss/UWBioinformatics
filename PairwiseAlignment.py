import Constants
from Logger import Logger
from Sequence import *

class PairwiseAlignment:
    """
    Computes a global/local pairwise alignment using dynamic programming algorithms

    Attributes:
        seqA    -- Sequence A
        seqB    -- Sequence B
        score   -- Scoring method
        path    -- Saves retrace path from dynamic algorithm

    Methods:
        Global  -- Needleman-Wunsch algorithm
        Local   -- Smith-Waterman algorithm
    """
    def __init__(self, seqA: Sequence, seqB: Sequence, score: Score):
        self.seqA = seqA.sequence
        self.seqB = seqB.sequence
        self.score = score
        self.path = []
        return

    def Global(self):
        # TODO
        
        return

    def Local(self):
        # TODO

        return
