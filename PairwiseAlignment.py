from Constants import *
from Logger import Logger
from Sequence import *
from Score import Score
from copy import deepcopy
import numpy as np

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
        self.alignments = None
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
            "h_gap": np.zeros((vlen, hlen)),
            "v_gap": np.zeros((vlen, hlen)),
            "match": np.zeros((vlen, hlen))
        }
        self.dp_array = np.zeros((vlen, hlen))
        self.direction = [["" for i in range(hlen)] for j in range(vlen)]

        # Initialize recurrence matrices
        matrices["v_gap"][0][0] = -np.Inf
        for i in range(1, hlen):
            matrices["v_gap"][0][i] = exist + i*extend
            self.dp_array[0][i] = exist + i*extend
            self.direction[0][i] += "L"

        matrices["h_gap"][0][0] = -np.Inf
        for j in range(1, vlen):
            matrices["h_gap"][j][0] = exist + j*extend
            self.dp_array[j][0] = exist + j*extend
            self.direction[j][0] += "U"

        matrices["match"][0][0] = -np.Inf

        # Affine gap penalty recurrence relation
        for j in range(1, vlen):
            for i in range(1, hlen):
                matrices["h_gap"][j][i] = max(matrices["h_gap"][j][i-1] + extend,
                                         self.dp_array[j][i-1] + exist + extend)
                matrices["v_gap"][j][i] = max(matrices["v_gap"][j-1][i] + extend,
                                         self.dp_array[j-1][i] + exist + extend)
                matrices["match"][j][i] = (self.dp_array[j-1][i-1]
                                           + matrix[self.vseq[j-1]][self.hseq[i-1]])
                self.dp_array[j][i] = max(matrices["h_gap"][j][i],
                                          matrices["v_gap"][j][i],
                                          matrices["match"][j][i])
                # Direction
                if self.dp_array[j][i] == matrices["h_gap"][j][i]:
                    self.direction[j][i] += "L"
                if self.dp_array[j][i] == matrices["v_gap"][j][i]:
                    self.direction[j][i] += "U"
                if self.dp_array[j][i] == matrices["match"][j][i]:
                    self.direction[j][i] += "D"

        # Find alignment(s)
        self.direction = np.array(self.direction)
        self.alignments = []
        self.__GlobalPaths("", hlen-1, vlen-1)
        for i in range(len(self.alignments)):
            print(f"Alignment {i+1}\n{self.alignments[i]}")
        print(f"Score: {self.dp_array[vlen-1][hlen-1]}")
        return

    def __GlobalPaths(self, temp_path: str, i: int, j: int):
        """This function was rough to code...

        Args:
            temp_path (list): Temporary path
            i (int): Optimal horizontal index in dp_array
            j (int): Optimal vertical index in dp_array
        """
        dir_dict = {
            "L": [-1, 0], # Left
            "U": [0, -1], # Up
            "D": [-1, -1] # Diagonal
        }

        # Perform depth first walk to origin while obtaining each matching
        if i > 0 or j > 0:
            curr = self.direction[j][i]
            for direction in curr:
                self.__GlobalPaths(temp_path + direction,
                            i + dir_dict[direction][0],
                            j + dir_dict[direction][1])
        else:
            temp_path = temp_path[::-1]
            hmatch = ""
            vmatch = ""
            hcounter = 0
            vcounter = 0
            for direction in temp_path:
                if direction == "D":
                    hmatch += self.hseq[hcounter]
                    vmatch += self.vseq[vcounter]
                    hcounter += 1
                    vcounter += 1
                if direction == "L":
                    hmatch += self.hseq[hcounter]
                    vmatch += "-"
                    hcounter += 1
                if direction == "U":
                    hmatch += "-"
                    vmatch += self.vseq[vcounter]
                    vcounter += 1
            match = hmatch + "\n" + vmatch
            self.alignments.append(match)
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
            "h_gap": np.zeros((vlen, hlen)),
            "v_gap": np.zeros((vlen, hlen)),
            "match": np.zeros((vlen, hlen))
        }
        self.dp_array = np.zeros((vlen, hlen))
        self.direction = [["" for i in range(hlen)] for j in range(vlen)]

        # Affine gap penalty recurrence relation
        for j in range(1, vlen):
            for i in range(1, hlen):
                matrices["h_gap"][j][i] = max(matrices["h_gap"][j][i-1] + extend,
                                         self.dp_array[j][i-1] + exist + extend, 0)
                matrices["v_gap"][j][i] = max(matrices["v_gap"][j-1][i] + extend,
                                         self.dp_array[j-1][i] + exist + extend, 0)
                matrices["match"][j][i] = (self.dp_array[j-1][i-1]
                                           + matrix[self.vseq[j-1]][self.hseq[i-1]])
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
        self.direction = np.array(self.direction)
        self.alignments = []
        optimum = np.max(self.dp_array)
        for j in range(vlen):
            for i in range(hlen):
                if self.dp_array[j][i] == optimum:
                    self.__LocalPaths("", i, j)
        for i in range(len(self.alignments)):
            print(f"Alignment {i+1}\n{self.alignments[i]}")
        print(f"Score: {optimum}")
        return
    
    def __LocalPaths(self, temp_path: str, i: int, j: int):
        """This one was easier...

        Args:
            temp_path (list): Temporary path
            i (int): Optimal horizontal index in dp_array
            j (int): Optimal vertical index in dp_array
        """
        dir_dict = {
            "L": [-1, 0], # Left
            "U": [0, -1], # Up
            "D": [-1, -1] # Diagonal
        }

        # Perform depth first walk to nearest 0 while obtaining matching
        if self.dp_array[j][i] != 0:
            curr = self.direction[j][i]
            for direction in curr:
                self.__LocalPaths(temp_path + direction,
                            i + dir_dict[direction][0],
                            j + dir_dict[direction][1])
        else:
            temp_path = temp_path[::-1]
            hmatch = ""
            vmatch = ""
            hcounter = 0
            vcounter = 0
            for direction in temp_path:
                if direction == "D":
                    hmatch += self.hseq[hcounter]
                    vmatch += self.vseq[vcounter]
                    hcounter += 1
                    vcounter += 1
                if direction == "L":
                    hmatch += self.hseq[hcounter]
                    vmatch += "-"
                    hcounter += 1
                if direction == "U":
                    hmatch += "-"
                    vmatch += self.vseq[vcounter]
                    vcounter += 1
            match = hmatch + "\n" + vmatch
            self.alignments.append(match)
        return

    def Summary(self):
        """Prints summary of pairwise alignment
        """
        print(self.dp_array)
        print(self.direction)
        return

score = Score(MatrixConstants.EXAMPLE_MATRIX, BlockConstants.NUCLEOTIDES, 0, -1)
print(score)
DNASeq1 = DNASequence("GTCGACGCA")
DNASeq2 = DNASequence("GATTACA")
alignment = PairwiseAlignment(DNASeq1, DNASeq2, score)
alignment.Local()
