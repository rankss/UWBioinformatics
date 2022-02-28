import numpy as np
from Sequence import Sequence
from Score import Score, BLOSUM62

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
    __DIR_DICT = {
        "L": [-1, 0], # Left
        "U": [0, -1], # Up
        "D": [-1, -1] # Diagonal
    }
    
    def __init__(self, sequence1: str, sequence2: str, score: Score):
        # self.logger = Logger()
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

    def __Path(self, temp_path: str, nodes: list, i: int, j: int, alignment=GLOBAL):
        if alignment == PairwiseAlignment.GLOBAL:
            comparison = i > 0 or j > 0
        if alignment == PairwiseAlignment.LOCAL:
            comparison = self.dpArray[j][i] != 0
            
        if comparison:
            curr = self.direction[j][i]
            for direction in curr:
                self.__Path(temp_path + direction, nodes + [(i, j)],
                            i + PairwiseAlignment.__DIR_DICT[direction][0],
                            j + PairwiseAlignment.__DIR_DICT[direction][1],
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
        hGap = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))
        vGap = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))
        match = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))

        # Initialize recurrence matrices
        vGap[0][0] = -np.Inf
        hGap[0][0] = -np.Inf
        match[0][0] = -np.Inf
        
        for i in range(1, hLen):
            vGap[0][i] = exist + i*extend
            self.dpArray[0][i] = exist + i*extend
            self.direction[0][i] += "L"

        for j in range(1, vLen):
            hGap[j][0] = exist + j*extend
            self.dpArray[j][0] = exist + j*extend
            self.direction[j][0] += "U"

        # Affine gap penalty recurrence relation
        for j in range(1, vLen):
            for i in range(1, hLen):
                hGap[j][i] = max(hGap[j][i-1] + extend,
                                        self.dpArray[j][i-1] + exist + extend)
                
                vGap[j][i] = max(vGap[j-1][i] + extend,
                                        self.dpArray[j-1][i] + exist + extend)
                
                match[j][i] = self.dpArray[j-1][i-1] + matrix[self.vSeq[j-1]][self.hSeq[i-1]]
                
                self.dpArray[j][i] = max(hGap[j][i], vGap[j][i], match[j][i])
                
                # Direction
                if self.dpArray[j][i] == hGap[j][i]:
                    self.direction[j][i] += "L"
                if self.dpArray[j][i] == vGap[j][i]:
                    self.direction[j][i] += "U"
                if self.dpArray[j][i] == match[j][i]:
                    self.direction[j][i] += "D"

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
        for j in range(1, vLen):
            for i in range(1, hLen):
                hGap[j][i] = max(0, hGap[j][i-1] + extend,
                                        self.dpArray[j][i-1] + exist + extend)
                
                vGap[j][i] = max(0, vGap[j-1][i] + extend,
                                        self.dpArray[j-1][i] + exist + extend)
                
                match[j][i] = self.dpArray[j-1][i-1] + matrix[self.vSeq[j-1]][self.hSeq[i-1]]
                
                self.dpArray[j][i] = max(0, hGap[j][i], vGap[j][i], match[j][i])
                
                # Direction
                if self.dpArray[j][i] == hGap[j][i]:
                    self.direction[j][i] += "L"
                if self.dpArray[j][i] == vGap[j][i]:
                    self.direction[j][i] += "U"
                if self.dpArray[j][i] == match[j][i]:
                    self.direction[j][i] += "D"

        # Create alignment(s)
        self.alignments = []
        self.paths = []
        self.optimal = np.max(self.dpArray)
        for j in range(vLen):
            for i in range(hLen):
                if self.dpArray[j][i] == self.optimal:
                    self.__Path("", [], i, j, PairwiseAlignment.LOCAL)
        return
    
    def Align(self, alignment=GLOBAL):
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
        __vGap = np.zeros((vLen, hLen + 1))
        __hGap = np.zeros((vLen, hLen + 1))
        __match = np.zeros((vLen, hLen + 1))

        # Initialize recurrence matrices
        __vGap[0][0] = -np.Inf
        __hGap[0][0] = -np.Inf
        __match[0][0] = -np.Inf
        
        for i in range(1, hLen):
            __vGap[0][i] = exist + i*extend
            dpArray[0][i] = exist + i*extend

        for j in range(1, vLen):
            __hGap[j][0] = exist + j*extend
            dpArray[j][0] = exist + j*extend

        # Affine gap penalty recurrence relation
        for j in range(1, vLen):
            for i in range(1, hLen):
                __hGap[j][i] = max(__hGap[j][i-1] + extend,
                                   dpArray[j][i-1] + exist + extend)
                
                __vGap[j][i] = max(__vGap[j-1][i] + extend,
                                   dpArray[j-1][i] + exist + extend)
                
                __match[j][i] = dpArray[j-1][i-1] + matrix[vSeq[j-1]][hSeq[i-1]]
                
                dpArray[j][i] = max(__hGap[j][i], __vGap[j][i], __match[j][i])
                
        return dpArray[vLen-1][hLen-1]

    @staticmethod
    def Local(hSeq: str, vSeq: str, score: Score) -> float:
        hSeq, vSeq = Sequence(hSeq), Sequence(vSeq)
        hLen, vLen = len(hSeq) + 1, len(vSeq) + 1
        exist, extend = score.existence, score.extension
        matrix = score.matrix
        
        # Define recurrence matrices
        dpArray = np.zeros((vLen, hLen))
        __vGap = np.zeros((vLen, hLen + 1))
        __hGap = np.zeros((vLen, hLen + 1))
        __match = np.zeros((vLen, hLen + 1))
        
        # Affine gap penalty recurrence relation
        for j in range(1, vLen):
            for i in range(1, hLen):
                __hGap[j][i] = max(0, __hGap[j][i-1] + extend,
                                   dpArray[j][i-1] + exist + extend)
                
                __vGap[j][i] = max(0, __vGap[j-1][i] + extend,
                                   dpArray[j-1][i] + exist + extend)
                
                __match[j][i] = dpArray[j-1][i-1] + matrix[vSeq[j-1]][hSeq[i-1]]
                
                dpArray[j][i] = max(0, __hGap[j][i], __vGap[j][i], __match[j][i])
                
        return np.max(dpArray)
        
    
class MultipleSequenceAlignment:
    """[summary]
    """
    def __init__(self, sequences: list, score: Score=BLOSUM62):
        self.sequences = sequences
        self.score = score
        self.count = len(self.sequences)
        return
    
    def ClustalW(self):
        """
          A  B  C  D  E
        A -  -  -  -  -
        B 9  -  -  -  -
        C 8  11 -  -  - 
        D 12 15 10 -  -
        E 15 18 13 5  -
        """
        # Create guide tree based on NW score
        guideMatrix = np.empty((self.count, self.count))
        for i, hSeq in enumerate(self.sequences[:-1]):
            for j, vSeq in enumerate(self.sequences, i + 1):
                guideMatrix[j][i] = PairwiseAlignment.Global(hSeq, vSeq, self.score)
        
        # UPGMA matrix merging
        symbols = {i:ord(i+65) for i in range(self.count)} # Define symbols
        maxIndex = np.unravel_index(guideMatrix.argmax(), guideMatrix.shape) # Find max index
        # Merge index and update matrix
        def updateMatrix(matrix: np.ndarray, mergeSymbols: tuple, iteration: int):
            guideMatrix = np.empty((self.count-1-iteration, self.count-1-iteration))
            
            
            pass
        
        return