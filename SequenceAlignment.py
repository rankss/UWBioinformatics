import numpy as np
from Sequence import Sequence
from Score import Score
from Error import InvalidAlignmentTypeError
from typing import Literal

class PairwiseAlignment:
    """Global/Local Pairwise Alignment

    Attributes:
        hseq (Sequence): Horizontal sequence
        vseq (Sequence): Vertical sequence
        score (Score): Scoring matrix with gap penalties
        dpArray (List): Contains final scoring of alignment
        direction (List): Contains direction per cell
        alignment (List): Contains all optimal alignments

    Methods:
        Global (None): Performs Needleman-Wunsch with affine gap penalty
        Local (None): Performs Smith-Waterman with affine gap penalty
        Summary (None): Returns summary of previous alignment
    """
    
    GLOBAL = 0
    LOCAL = 1
    __VALID_ALIGNMENT_TYPE = {GLOBAL, LOCAL}
    __DIR_DICT = {
        "L": [-1, 0], # Left
        "U": [0, -1], # Up
        "D": [-1, -1] # Diagonal
    }
    
    def __init__(self, sequence1: str, sequence2: str, score: Score):
        self.score = score
        self.hSeq = Sequence(sequence1)
        self.vSeq = Sequence(sequence2)
        self.dpArray = None
        self.direction = None
        self.paths = None
        self.alignments = None
        self.optimal = None

        return

    def __Match(self, temp_path: str): # Match has bug for local alignment
        hMatch, vMatch = "", ""
        hCounter, vCounter = 0, 0
        for direction in temp_path:
            if direction == "D":
                hMatch += self.hSeq[hCounter]
                vMatch += self.vSeq[vCounter]
                hCounter += 1
                vCounter += 1
            if direction == "L":
                hMatch += self.hSeq[hCounter]
                vMatch += "-"
                hCounter += 1
            if direction == "U":
                hMatch += "-"
                vMatch += self.vSeq[vCounter]
                vCounter += 1
        match = f"{hMatch}\n{vMatch}"
        self.alignments.append(match)
        return

    def __Path(self, temp_path: str, nodes: list, col: int, row: int, alignment: Literal=GLOBAL):
        if alignment == PairwiseAlignment.GLOBAL:
            comparison = col > 0 or row > 0
        if alignment == PairwiseAlignment.LOCAL:
            comparison = self.dpArray[row, col] != 0
            
        if comparison:
            curr = self.direction[row][col]
            for direction in curr:
                self.__Path(temp_path + direction, nodes + [(col, row)],
                            col + PairwiseAlignment.__DIR_DICT[direction][0],
                            row + PairwiseAlignment.__DIR_DICT[direction][1],
                            alignment)
        else:
            self.__Match(temp_path[::-1])
            self.paths.append(nodes)
        return
        
    def __Global(self):
        """Performs Needleman-Wunsch for sequence pair with affine gap penalty
        """
        hLen, vLen = len(self.hSeq) + 1, len(self.vSeq) + 1
        exist, extend = self.score.existence, self.score.extension
        matrix = self.score.matrix
        
        # Define recurrence matrices
        self.dpArray = np.zeros((vLen, hLen))
        self.direction = [["" for i in range(hLen)] for j in range(vLen)]
        hGap = np.full((len(self.vSeq) + 1, len(self.hSeq) + 1), -np.inf)
        vGap = np.full((len(self.vSeq) + 1, len(self.hSeq) + 1), -np.inf)
        match = np.full((len(self.vSeq) + 1, len(self.hSeq) + 1), -np.inf)

        
        for col in range(1, hLen):
            vGap[0, col] = exist + col*extend
            self.dpArray[0, col] = exist + col*extend
            self.direction[0][col] += "L"

        for row in range(1, vLen):
            hGap[row, 0] = exist + row*extend
            self.dpArray[row, 0] = exist + row*extend
            self.direction[row][0] += "U"

        # Affine gap penalty recurrence relation
        for row in range(1, vLen):
            for col in range(1, hLen):
                hGap[row, col] = max(hGap[row, col-1] + extend,
                                     self.dpArray[row, col-1] + exist + extend)
                
                vGap[row, col] = max(vGap[row-1, col] + extend,
                                     self.dpArray[row-1, col] + exist + extend)
                
                match[row, col] = self.dpArray[row-1, col-1] + matrix[self.vSeq[row-1]][self.hSeq[col-1]]
                
                self.dpArray[row, col] = max(hGap[row, col], vGap[row, col], match[row, col])
                
                # Direction
                if self.dpArray[row, col] == hGap[row, col]:
                    self.direction[row][col] += "L"
                if self.dpArray[row, col] == vGap[row, col]:
                    self.direction[row][col] += "U"
                if self.dpArray[row, col] == match[row, col]:
                    self.direction[row][col] += "D"

        # Find alignment(s)
        self.optimal = self.dpArray[vLen-1][hLen-1]
        self.alignments = []
        self.paths = []
        self.__Path("", [], hLen-1, vLen-1, PairwiseAlignment.GLOBAL)
        return 

    def __Local(self):
        """Performs Smith-Waterman for sequence pair with affine gap penalty
        """
        hLen, vLen = len(self.hSeq) + 1, len(self.vSeq) + 1
        exist, extend = self.score.existence, self.score.extension
        matrix = self.score.matrix

        # Define recurrence matrices
        self.dpArray = np.zeros((vLen, hLen))
        self.direction = [["" for i in range(hLen)] for j in range(vLen)]
        hGap = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))
        vGap = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))
        match = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))

        # Affine gap penalty recurrence relation
        for row in range(1, vLen):
            for col in range(1, hLen):
                hGap[row, col] = max(0, hGap[row, col-1] + extend,
                                     self.dpArray[row, col-1] + exist + extend)
                
                vGap[row, col] = max(0, vGap[row-1, col] + extend,
                                     self.dpArray[row-1, col] + exist + extend)
                
                match[row, col] = self.dpArray[row-1, col-1] + matrix[self.vSeq[row-1]][self.hSeq[col-1]]
                
                self.dpArray[row, col] = max(0, hGap[row, col], vGap[row, col], match[row, col])
                
                # Direction
                if self.dpArray[row, col] == hGap[row, col]:
                    self.direction[row][col] += "L"
                if self.dpArray[row, col] == vGap[row, col]:
                    self.direction[row][col] += "U"
                if self.dpArray[row, col] == match[row, col]:
                    self.direction[row][col] += "D"

        # Create alignment(s)
        self.alignments = []
        self.paths = []
        self.optimal = np.max(self.dpArray)
        for row in range(vLen):
            for col in range(hLen):
                if self.dpArray[row, col] == self.optimal:
                    self.__Path("", [], col, row, PairwiseAlignment.LOCAL)
        return
    
    def Align(self, alignment: Literal=GLOBAL):
        if alignment not in PairwiseAlignment.__VALID_ALIGNMENT_TYPE:
            raise InvalidAlignmentTypeError()
    
        if alignment == PairwiseAlignment.GLOBAL:
            self.__Global()
        if alignment == PairwiseAlignment.LOCAL:
            self.__Local()
        return
    
    def Info(self):
        pass
    
    @staticmethod
    def Global(hSeq: str, vSeq: str, score: Score) -> float:
        hSeq, vSeq = Sequence(hSeq), Sequence(vSeq)
        hLen, vLen = len(hSeq) + 1, len(vSeq) + 1
        exist, extend = score.existence, score.extension
        matrix = score.matrix
        
        # Define recurrence matrices
        dpArray = np.zeros((vLen, hLen))
        vGap = np.full((vLen, hLen + 1), -np.inf)
        hGap = np.full((vLen, hLen + 1), -np.inf)
        match = np.full((vLen, hLen + 1), -np.inf)
        
        for col in range(1, hLen):
            vGap[0, col] = exist + col*extend
            dpArray[0, col] = exist + col*extend

        for row in range(1, vLen):
            hGap[row, 0] = exist + row*extend
            dpArray[row, 0] = exist + row*extend

        # Affine gap penalty recurrence relation
        for row in range(1, vLen):
            for col in range(1, hLen):
                hGap[row, col] = max(hGap[row, col-1] + extend,
                                     dpArray[row, col-1] + exist + extend)
                
                vGap[row, col] = max(vGap[row-1, col] + extend,
                                     dpArray[row-1, col] + exist + extend)
                
                match[row, col] = dpArray[row-1, col-1] + matrix[vSeq[row-1]][hSeq[col-1]]
                
                dpArray[row, col] = max(hGap[row, col], vGap[row, col], match[row, col])
                
        return dpArray[vLen-1][hLen-1]

    @staticmethod
    def Local(hSeq: str, vSeq: str, score: Score) -> float:
        hSeq, vSeq = Sequence(hSeq), Sequence(vSeq)
        hLen, vLen = len(hSeq) + 1, len(vSeq) + 1
        exist, extend = score.existence, score.extension
        matrix = score.matrix
        
        # Define recurrence matrices
        dpArray = np.zeros((vLen, hLen))
        vGap = np.zeros((vLen, hLen + 1))
        hGap = np.zeros((vLen, hLen + 1))
        match = np.zeros((vLen, hLen + 1))
        
        # Affine gap penalty recurrence relation
        for row in range(1, vLen):
            for col in range(1, hLen):
                hGap[row, col] = max(0, hGap[row, col-1] + extend,
                                 dpArray[row, col-1] + exist + extend)
                
                vGap[row, col] = max(0, vGap[row-1, col] + extend,
                                 dpArray[row-1, col] + exist + extend)
                
                match[row, col] = dpArray[row-1, col-1] + matrix[vSeq[row-1]][hSeq[col-1]]
                
                dpArray[row, col] = max(0, hGap[row, col], vGap[row, col], match[row, col])
                
        return np.max(dpArray)
        
    
class MultipleSequenceAlignment:
    """[summary]
    """
    def __init__(self, sequences: list, score: Score):
        self.sequences = sequences
        self.score = score
        self.count = len(self.sequences)
        return
    
    def ClustalW(self):
        return