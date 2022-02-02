from Constants import *
from Logger import Logger
from Sequence import *
from Score import Score
from math import inf
from copy import deepcopy

class PairwiseAlignment:
    """Global/Local Pairwise Alignment
    
    Attributes:
        hseq (Sequence): Horizontal sequence
        vseq (Sequence): Vertical sequence
        score (Score): Scoring matrix with gap penalties
        dp_array (List): Contains final scoring of alignment
        direction (List): Contains direction per cell
        alignment (List): Contains all optimal alignments
        
    Methods:
        Global (None): Performs Needleman-Wunsch with affine gap penalty
        Local (None): Performs Smith-Waterman with affine gap penalty
        Summary (None): Returns summary of previous alignment
    """
    def __init__(self, sequence1: Sequence, sequence2: Sequence, score: Score):
        # self.logger = Logger()
        self.hseq = sequence1.sequence
        self.vseq = sequence2.sequence
        self.score = score
        self.dp_array = None
        self.direction = None
        self.alignment = None
        return

    def Global(self):
        """Performs Needleman-Wunsch for sequence pair with affine gap penalty
        """
        hlen = len(self.hseq) + 1
        vlen = len(self.vseq) + 1
        exist = self.score.existence
        extend = self.score.extension
        matrix = self.score.matrix

        # Define recurrence matrices
        matrices = {
            "h_gap": [[0 for i in range(hlen)] for j in range(vlen)],
            "v_gap": [[0 for i in range(hlen)] for j in range(vlen)],
            "match": [[0 for i in range(hlen)] for j in range(vlen)]
        }
        self.dp_array = [[0 for i in range(hlen)] for j in range(vlen)]
        self.direction = [["" for i in range(hlen)] for j in range(vlen)]

        # Initialize recurrence matrices
        matrices["v_gap"][0][0] = -inf
        for i in range(1, hlen):
            matrices["v_gap"][0][i] = exist + i*extend
            self.dp_array[0][i] = exist + i*extend
            self.direction[0][i] += "L"

        matrices["h_gap"][0][0] = -inf
        for j in range(1, vlen):
            matrices["h_gap"][j][0] = exist + j*extend
            self.dp_array[j][0] = exist + j*extend
            self.direction[j][0] += "U"

        matrices["match"][0][0] = -inf

        # Affine gap penalty recurrence relation
        for j in range(1, vlen):
            for i in range(1, hlen):
                matrices["h_gap"][j][i] = max(matrices["h_gap"][j][i-1] + extend,
                                         self.dp_array[j][i-1] + exist + extend)
                matrices["v_gap"][j][i] = max(matrices["v_gap"][j-1][i] + extend,
                                         self.dp_array[j-1][i] + exist + extend)
                matrices["match"][j][i] = self.dp_array[j-1][i-1] + matrix[self.vseq[j-1]][self.hseq[i-1]]
                self.dp_array[j][i] = max(matrices["v_gap"][j][i],
                                          matrices["h_gap"][j][i],
                                          matrices["match"][j][i])
                # Direction
                if self.dp_array[j][i] == matrices["h_gap"][j][i]:
                    self.direction[j][i] += "L"
                if self.dp_array[j][i] == matrices["v_gap"][j][i]:
                    self.direction[j][i] += "U"
                if self.dp_array[j][i] == matrices["match"][j][i]:
                    self.direction[j][i] += "D"

        # Find alignment(s)
        # TODO
        
        return

    def Local(self):
        """Performs Smith-Waterman for sequence pair with affine gap penalty
        """
        hlen = len(self.hseq) + 1
        vlen = len(self.vseq) + 1
        exist = self.score.existence
        extend = self.score.extension
        matrix = self.score.matrix

        # Define recurrence matrices
        matrices = {
            "h_gap": [[0 for i in range(hlen)] for j in range(vlen)],
            "v_gap": [[0 for i in range(hlen)] for j in range(vlen)],
            "match": [[0 for i in range(hlen)] for j in range(vlen)]
        }
        self.dp_array = [[0 for i in range(hlen)] for j in range(vlen)]
        self.direction = [["" for i in range(hlen)] for j in range(vlen)]

        # Affine gap penalty recurrence relation
        for j in range(1, vlen):
            for i in range(1, hlen):
                matrices["h_gap"][j][i] = max(matrices["h_gap"][j][i-1] + extend,
                                         self.dp_array[j][i-1] + exist + extend, 0)
                matrices["v_gap"][j][i] = max(matrices["v_gap"][j-1][i] + extend,
                                         self.dp_array[j-1][i] + exist + extend, 0)
                matrices["match"][j][i] = self.dp_array[j-1][i-1] + matrix[self.vseq[j-1]][self.hseq[i-1]]
                self.dp_array[j][i] = max(matrices["v_gap"][j][i],
                                          matrices["h_gap"][j][i],
                                          matrices["match"][j][i], 0)
                # Direction
                if self.dp_array[j][i] == matrices["h_gap"][j][i]:
                    self.direction[j][i] += "L"
                if self.dp_array[j][i] == matrices["v_gap"][j][i]:
                    self.direction[j][i] += "U"
                if self.dp_array[j][i] == matrices["match"][j][i]:
                    self.direction[j][i] += "D"

        # Create alignment(s)
        # TODO
        
        return

    def Summary(self):
        """Prints summary of pairwise alignment
        """
        print(self.score)
        # TODO

score = Score(MatrixConstants.EXAMPLE_MATRIX, BlockConstants.NUCLEOTIDES, 0, 2)
print(score)
DNASeq1 = DNASequence("CAAGAC")
DNASeq2 = DNASequence("GAAC")
alignment = PairwiseAlignment(DNASeq1, DNASeq2, score)
alignment.Global()