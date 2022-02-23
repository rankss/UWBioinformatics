import numpy as np
import matplotlib.pyplot as plt
from Sequence import Sequence
from Matrix import EXAMPLE_MATRIX, Score, BLOSUM62, EXAMPLE_MATRIX

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
    def __init__(self, sequence1: str, sequence2: str, score: Score=EXAMPLE_MATRIX):
        # self.logger = Logger()
        self.score = score
        self.hSeq = Sequence(sequence1, self.score.sequenceType)
        self.vSeq = Sequence(sequence2, self.score.sequenceType)
        self.dpArray = None
        self.direction = None
        self.paths = None
        self.alignments = None
        self.optimal = None
        
        self.__hGap = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))
        self.__vGap = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))
        self.__match = np.zeros((len(self.vSeq) + 1, len(self.hSeq) + 1))
        
        self.__dir_dict = {
            "L": [-1, 0], # Left
            "U": [0, -1], # Up
            "D": [-1, -1] # Diagonal
        }
        self.__colours = {
            1: "#ADD8E6", # Blue
            2: "#90EE90", # Green
            3: "#FED8B1", # Orange
            4: "#F4C2C2", # Pink
            5: "#FFFCBB", # C3
            6: "#B19CD9", # C4
            7: "#F7D8BA"  # C5
        }

        return
    
    def OptimalScore(self):
        return self.optimal

    def __Match(self, temp_path: str):
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

    def Global(self):
        """Performs Needleman-Wunsch for sequence pair with affine gap penalty
        """
        hLen, vLen = len(self.hSeq) + 1, len(self.vSeq) + 1
        exist, extend = self.score.existence, self.score.extension
        matrix = self.score.matrix
        
        # Define recurrence matrices
        self.dpArray = np.zeros((vLen, hLen))
        self.direction = [["" for i in range(hLen)] for j in range(vLen)]

        # Initialize recurrence matrices
        self.__vGap[0][0] = -np.Inf
        self.__hGap[0][0] = -np.Inf
        self.__match[0][0] = -np.Inf
        
        for i in range(1, hLen):
            self.__vGap[0][i] = exist + i*extend
            self.dpArray[0][i] = exist + i*extend
            self.direction[0][i] += "L"

        for j in range(1, vLen):
            self.__hGap[j][0] = exist + j*extend
            self.dpArray[j][0] = exist + j*extend
            self.direction[j][0] += "U"

        # Affine gap penalty recurrence relation
        for j in range(1, vLen):
            for i in range(1, hLen):
                self.__hGap[j][i] = max(self.__hGap[j][i-1] + extend,
                                        self.dpArray[j][i-1] + exist + extend)
                
                self.__vGap[j][i] = max(self.__vGap[j-1][i] + extend,
                                        self.dpArray[j-1][i] + exist + extend)
                
                self.__match[j][i] = self.dpArray[j-1][i-1] + matrix[self.vSeq[j-1]][self.hSeq[i-1]]
                
                self.dpArray[j][i] = max(self.__hGap[j][i], self.__vGap[j][i], self.__match[j][i])
                
                # Direction
                if self.dpArray[j][i] == self.__hGap[j][i]:
                    self.direction[j][i] += "L"
                if self.dpArray[j][i] == self.__vGap[j][i]:
                    self.direction[j][i] += "U"
                if self.dpArray[j][i] == self.__match[j][i]:
                    self.direction[j][i] += "D"
                    
        def GlobalPaths(temp_path: str, nodes: list, i: int, j: int):
            if i > 0 or j > 0:
                curr = self.direction[j][i]
                for direction in curr:
                    GlobalPaths(temp_path + direction, nodes + [[i, j]],
                                i + self.__dir_dict[direction][0],
                                j + self.__dir_dict[direction][1])
            else:
                self.__Match(temp_path[::-1])
                self.paths.append(nodes)
            return

        # Find alignment(s)
        self.optimal = self.dpArray[vLen-1][hLen-1]
        self.alignments = []
        self.paths = []
        GlobalPaths("", [], hLen-1, vLen-1)
        return

    def Local(self):
        """Performs Smith-Waterman for sequence pair with affine gap penalty
        """
        hLen, vLen = len(self.hSeq) + 1, len(self.vSeq) + 1
        exist, extend = self.score.existence, self.score.extension
        matrix = self.score.matrix

        # Define recurrence matrices
        self.dpArray = np.zeros((vLen, hLen))
        self.direction = [["" for i in range(hLen)] for j in range(vLen)]

        # Affine gap penalty recurrence relation
        for j in range(1, vLen):
            for i in range(1, hLen):
                self.__hGap[j][i] = max(0, self.__hGap[j][i-1] + extend,
                                        self.dpArray[j][i-1] + exist + extend)
                
                self.__vGap[j][i] = max(0, self.__vGap[j-1][i] + extend,
                                        self.dpArray[j-1][i] + exist + extend)
                
                self.__match[j][i] = self.dpArray[j-1][i-1] + matrix[self.vSeq[j-1]][self.hSeq[i-1]]
                
                self.dpArray[j][i] = max(0, self.__hGap[j][i], self.__vGap[j][i], self.__match[j][i])
                
                # Direction
                if self.dpArray[j][i] == self.__hGap[j][i]:
                    self.direction[j][i] += "L"
                if self.dpArray[j][i] == self.__vGap[j][i]:
                    self.direction[j][i] += "U"
                if self.dpArray[j][i] == self.__match[j][i]:
                    self.direction[j][i] += "D"
                    
        def LocalPaths(temp_path: str, nodes: list, i: int, j: int):
            # Perform depth first walk to nearest 0 while obtaining matching
            if self.dpArray[j][i] != 0:
                curr = self.direction[j][i]
                for direction in curr:
                    LocalPaths(temp_path + direction, nodes + [[i, j]],
                               i + self.__dir_dict[direction][0],
                               j + self.__dir_dict[direction][1])
            else:
                self.__Match(temp_path[::-1])
                self.paths.append(nodes)
            return

        # Create alignment(s)
        self.alignments = []
        self.paths = []
        self.optimal = np.max(self.dpArray)
        for j in range(vLen):
            for i in range(hLen):
                if self.dpArray[j][i] == self.optimal:
                    LocalPaths("", [], i, j)
        print(self.paths)
        for i in range(len(self.alignments)):
            print(f"Alignment {i+1}\n{self.alignments[i]}")
        print(f"Score: {self.optimal}")
        return
    
    def Summary(self):
        """Prints summary of pairwise alignment
        """
        row_labels = [" "] + list(self.vSeq.sequence)
        col_labels = [" "] + list(self.hSeq.sequence)

        cell_colours = [["w" for i in range(len(col_labels))] for j in range(len(row_labels))]

        for i in range(len(self.paths)):
            for j in range(len(self.paths[i])):
                node = self.paths[i][j]
                if j == 0:
                    cell_colours[node[1]][node[0]] = self.__colours[1]
                else:
                    if cell_colours[node[1]][node[0]] == "w":
                        cell_colours[node[1]][node[0]] = self.__colours[i+3]
                    elif cell_colours[node[1]][node[0]] == self.__colours[1]:
                        pass
                    else:
                        cell_colours[node[1]][node[0]] = self.__colours[2]

        fig, axis = plt.subplots(2)
        fig.suptitle("Trace")
        axis[0].set_axis_off()
        axis[0].table(cellText = self.dpArray,
                      cellColours = cell_colours,
                      rowLabels = row_labels,
                      rowColours = ["#ADD8E6" for row in row_labels],
                      rowLoc = "center",
                      colLabels = col_labels,
                      colColours = ["#ADD8E6" for col in col_labels],
                      colWidths = [0.05 for col in col_labels],
                      colLoc = "center",
                      cellLoc = "center",
                      loc = "center")

        axis[1].set_axis_off()
        axis[1].table(cellText = self.direction,
                      cellColours = cell_colours,
                      rowLabels = row_labels,
                      rowColours = ["#ADD8E6" for row in row_labels],
                      rowLoc = "center",
                      colLabels = col_labels,
                      colColours = ["#ADD8E6" for col in col_labels],
                      colWidths = [0.05 for col in col_labels],
                      colLoc = "center",
                      cellLoc = "center",
                      loc = "center")
        plt.show()
        return

class MultipleSequenceAlignment:
    """[summary]
    """
    def __init__(self, sequences: list, sequenceType: list, Score: Score=BLOSUM62):
        self.sequences = sequences
        return
    
    def ClustalW(self):
        # Construct NW for all pairs
        # Create guide tree based on score
        for i, sequence1 in enumerate(self.sequences):
            for sequence2 in self.sequences[i:]:
                pass
        
        return