import numpy as np
import matplotlib.pyplot as plt
from Sequence import Sequence
from Matrix import Score

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
    def __init__(self, sequence1: str, sequence2: str, score: Score):
        # self.logger = Logger()
        self.score = score
        self.hseq = Sequence(sequence1, self.score.sequenceType)
        self.vseq = Sequence(sequence2, self.score.sequenceType)
        self.dpArray = None
        self.direction = None
        self.paths = None
        self.alignments = None

        self.__matrices = {
            "hgap": np.zeros((len(self.vseq.sequence)+1, len(self.hseq.sequence)+1)),
            "vgap": np.zeros((len(self.vseq.sequence)+1, len(self.hseq.sequence)+1)),
            "match": np.zeros((len(self.vseq.sequence)+1, len(self.hseq.sequence)+1))
        }
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

    def __Match(self, temp_path: str):
        hseq = self.hseq.sequence
        vseq = self.vseq.sequence
        hmatch = ""
        vmatch = ""
        hcounter = 0
        vcounter = 0
        for direction in temp_path:
            if direction == "D":
                hmatch += hseq[hcounter]
                vmatch += vseq[vcounter]
                hcounter += 1
                vcounter += 1
            if direction == "L":
                hmatch += hseq[hcounter]
                vmatch += "-"
                hcounter += 1
            if direction == "U":
                hmatch += "-"
                vmatch += vseq[vcounter]
                vcounter += 1
        match = hmatch + "\n" + vmatch
        self.alignments.append(match)
        return

    def Global(self):
        """Performs Needleman-Wunsch for sequence pair with affine gap penalty
        """
        hseq = self.hseq.sequence
        vseq = self.vseq.sequence
        hlen = len(hseq) + 1
        vlen = len(vseq) + 1
        exist = self.score.existence
        extend = self.score.extension
        matrix = self.score.matrix

        # Define recurrence matrices
        self.dpArray = np.zeros((vlen, hlen))
        self.direction = [["" for i in range(hlen)] for j in range(vlen)]

        # Initialize recurrence matrices
        self.__matrices["vgap"][0][0] = -np.Inf
        for i in range(1, hlen):
            self.__matrices["vgap"][0][i] = exist + i*extend
            self.dpArray[0][i] = exist + i*extend
            self.direction[0][i] += "L"

        self.__matrices["hgap"][0][0] = -np.Inf
        for j in range(1, vlen):
            self.__matrices["hgap"][j][0] = exist + j*extend
            self.dpArray[j][0] = exist + j*extend
            self.direction[j][0] += "U"

        self.__matrices["match"][0][0] = -np.Inf

        # Affine gap penalty recurrence relation
        for j in range(1, vlen):
            for i in range(1, hlen):
                self.__matrices["hgap"][j][i] = max(self.__matrices["hgap"][j][i-1] + extend,
                                                     self.dpArray[j][i-1] + exist + extend)
                self.__matrices["vgap"][j][i] = max(self.__matrices["vgap"][j-1][i] + extend,
                                                     self.dpArray[j-1][i] + exist + extend)
                self.__matrices["match"][j][i] = self.dpArray[j-1][i-1] + matrix[vseq[j-1]][hseq[i-1]]
                self.dpArray[j][i] = max(self.__matrices["hgap"][j][i],
                                          self.__matrices["vgap"][j][i],
                                          self.__matrices["match"][j][i])
                # Direction
                if self.dpArray[j][i] == self.__matrices["hgap"][j][i]:
                    self.direction[j][i] += "L"
                if self.dpArray[j][i] == self.__matrices["vgap"][j][i]:
                    self.direction[j][i] += "U"
                if self.dpArray[j][i] == self.__matrices["match"][j][i]:
                    self.direction[j][i] += "D"

        # Find alignment(s)
        self.direction = np.array(self.direction)
        self.alignments = []
        self.paths = []
        self.__GlobalPaths("", [], hlen-1, vlen-1)
        for i in range(len(self.alignments)):
            print(f"Alignment {i+1}\n{self.alignments[i]}")
        print(f"Score: {self.dpArray[vlen-1][hlen-1]}")
        return

    def __GlobalPaths(self, temp_path: str, nodes: list, i: int, j: int):
        """This function was rough to code...

        Args:
            temp_path (list): Temporary path
            nodes (list): Traversed nodes
            i (int): Optimal horizontal index in dpArray
            j (int): Optimal vertical index in dpArray
        """
        # Perform depth first walk to origin while obtaining each matching
        if i > 0 or j > 0:
            curr = self.direction[j][i]
            for direction in curr:
                self.__GlobalPaths(temp_path + direction, nodes + [[i, j]],
                                   i + self.__dir_dict[direction][0],
                                   j + self.__dir_dict[direction][1])
        else:
            self.__Match(temp_path[::-1])
            self.paths.append(nodes)
        return

    def Local(self):
        """Performs Smith-Waterman for sequence pair with affine gap penalty
        """
        hseq = self.hseq.sequence
        vseq = self.vseq.sequence
        hlen = len(hseq) + 1
        vlen = len(vseq) + 1
        exist = self.score.existence
        extend = self.score.extension
        matrix = self.score.matrix

        # Define recurrence matrices
        self.dpArray = np.zeros((vlen, hlen))
        self.direction = [["" for i in range(hlen)] for j in range(vlen)]

        # Affine gap penalty recurrence relation
        for j in range(1, vlen):
            for i in range(1, hlen):
                self.__matrices["hgap"][j][i] = max(self.__matrices["hgap"][j][i-1] + extend,
                                                    self.dpArray[j][i-1] + exist + extend, 0)
                self.__matrices["vgap"][j][i] = max(self.__matrices["vgap"][j-1][i] + extend,
                                                    self.dpArray[j-1][i] + exist + extend, 0)
                self.__matrices["match"][j][i] = self.dpArray[j-1][i-1] + matrix[vseq[j-1]][hseq[i-1]]
                self.dpArray[j][i] = max(self.__matrices["vgap"][j][i],
                                         self.__matrices["hgap"][j][i],
                                         self.__matrices["match"][j][i], 0)
                # Direction
                if self.dpArray[j][i] == self.__matrices["hgap"][j][i]:
                    self.direction[j][i] += "L"
                if self.dpArray[j][i] == self.__matrices["vgap"][j][i]:
                    self.direction[j][i] += "U"
                if self.dpArray[j][i] == self.__matrices["match"][j][i]:
                    self.direction[j][i] += "D"

        # Create alignment(s)
        self.direction = np.array(self.direction)
        self.alignments = []
        self.paths = []
        optimum = np.max(self.dpArray)
        for j in range(vlen):
            for i in range(hlen):
                if self.dpArray[j][i] == optimum:
                    self.__LocalPaths("", [], i, j)
        for i in range(len(self.alignments)):
            print(f"Alignment {i+1}\n{self.alignments[i]}")
        print(f"Score: {optimum}")
        return

    def __LocalPaths(self, temp_path: str, nodes: list, i: int, j: int):
        """This one was easier...

        Args:
            temp_path (list): Temporary path
            i (int): Optimal horizontal index in dpArray
            j (int): Optimal vertical index in dpArray
        """
        # Perform depth first walk to nearest 0 while obtaining matching
        if self.dpArray[j][i] != 0:
            curr = self.direction[j][i]
            for direction in curr:
                self.__LocalPaths(temp_path + direction, nodes + [[i, j]],
                                  i + self.__dir_dict[direction][0],
                                  j + self.__dir_dict[direction][1])
        else:
            self.__Match(temp_path[::-1])
            self.paths.append(nodes)
        return

    def Summary(self):
        """Prints summary of pairwise alignment
        """
        row_labels = [" "] + list(self.vseq.sequence)
        col_labels = [" "] + list(self.hseq.sequence)

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
    def __init__(self, sequences: list):
        self.sequences = sequences
        return
