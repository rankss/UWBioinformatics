from Constants import *
from Logger import Logger
from Sequence import *
from Score import Score
from math import inf

class PairwiseAlignment:
    """[summary]
    """
    def __init__(self, sequence1: Sequence, sequence2: Sequence, score: Score):
        self.logger = Logger()
        self.hsequence = sequence1.sequence
        self.vsequence = sequence2.sequence
        self.score = score
        self.dp_array = None
        self.direction = None
        self.alignment = None
        return

    def __str__(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        # TODO
        return ""

    def Global(self):
        """Performs Needleman-Wunsch for sequence pair with affine gap penalty
        """
        hseq = self.hsequence
        vseq = self.vsequence
        hlen = len(hseq)
        vlen = len(vseq)
        exist = self.score.existence
        extend = self.score.extension
        matrix = self.score.matrix

        # Define recurrence matrices
        h_gap_matrix = [[0 for i in range(hlen+1)] for j in range(vlen+1)]
        v_gap_matrix = [[0 for i in range(hlen+1)] for j in range(vlen+1)]
        match_matrix = [[0 for i in range(hlen+1)] for j in range(vlen+1)]
        self.dp_array = [[0 for i in range(hlen+1)] for j in range(vlen+1)]
        self.direction = [["" for i in range(hlen+1)] for j in range(vlen+1)]

        # Initialize recurrence matrices
        v_gap_matrix[0][0] = -inf
        for i in range(1, hlen+1):
            v_gap_matrix[0][i] = exist + i*extend
            self.dp_array[0][i] = exist + i*extend

        h_gap_matrix[0][0] = -inf
        for j in range(1, vlen+1):
            h_gap_matrix[j][0] = exist + j*extend
            self.dp_array[j][0] = exist + j*extend

        match_matrix[0][0] = -inf

        # Affine gap penalty recurrence relation
        for j in range(1, vlen+1):
            for i in range(1, hlen+1):
                h_gap_matrix[j][i] = max(h_gap_matrix[j][i-1] + extend,
                                         self.dp_array[j][i-1] + exist + extend)
                v_gap_matrix[j][i] = max(v_gap_matrix[j-1][i] + extend,
                                         self.dp_array[j-1][i] + exist + extend)
                match_matrix[j][i] = self.dp_array[j-1][i-1] + matrix[vseq[j-1]][hseq[i-1]]
                self.dp_array[j][i] = max(v_gap_matrix[j][i], 
                                          h_gap_matrix[j][i],
                                          match_matrix[j][i])
                # Direction
                if self.dp_array[j][i] == h_gap_matrix[j][i]:
                    self.direction[j][i] += "L"
                if self.dp_array[j][i] == v_gap_matrix[j][i]:
                    self.direction[j][i] += "U"
                if self.dp_array[j][i] == match_matrix[j][i]:
                    self.direction[j][i] += "D"

        return

    def Local(self):
        """Performs Smith-Waterman for sequence pair
        """
        # TODO

        return

    def __Trace(self, dp_array: list, i: int, j: int):
        """[summary]

        Args:
            dp_array (list): [description]
            i (int): [description]
            j (int): [description]
        """
        # TODO
        
        return

score = Score(MatrixConstants.EXAMPLE_MATRIX, BlockConstants.NUCLEOTIDES, 0, 2)
DNASeq1 = DNASequence("CAAGAC")
DNASeq2 = DNASequence("GAAC")
alignment = PairwiseAlignment(DNASeq1, DNASeq2, score)
alignment.Global()