from Sequence import *

class Score:
    """
    Scoring scheme for pairwise alignment

    Attributes:
        matrix      -- Scoring matrix
        existence   -- Value for gap existence
        extension   -- Value for gap extension

    Methods:
        Matrix      -- Prints scoring matrix according to format indicated
    """

    def __init__(self, matrix: dict, existence: int, extension: int):
        self.matrix = matrix
        self.existence = existence
        self.extension = extension
        return

    def Matrix(self):
        """ 
        Output Format
        =====================
        |   | A | G | C | T |
        |---|---|---|---|---|
        | A | 1 |-1 |-1 |-1 |
        |---|---|---|---|---|
        | G |-1 | 1 |-1 |-1 |
        |---|---|---|---|---|
        | C |-1 |-1 | 1 |-1 |
        |---|---|---|---|---|
        | T |-1 |-1 |-1 | 1 |
        =====================
        """
        sequence = ['A', 'G', 'C', 'T']
        print("=====================")
        print("|   | A | G | C | T |")
        for vnt in sequence:
            values = f"|{vnt.rjust(2).ljust(3)}|"
            for hnt in sequence:
                values += f"{str(self.matrix[vnt][hnt]).rjust(2).ljust(3)}|"
            print("|---|---|---|---|---|")
            print(values)
        print("=====================")
        return

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
