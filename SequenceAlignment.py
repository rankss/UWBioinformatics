import numpy as np
import matplotlib.pyplot as plt
from Constants import Block
from Error import InvalidSequenceError, InvalidMatrixError

__all__ = ["Score", "PairwiseAlignment"]

class Sequence:
    """Base class for sequence strings

    Attributes:
        sequence    -- Stores sequence information after stripping and converting to upper case

    Methods:
        Validate    -- Implemented by subclasses
        Clean       -- Basic cleaning of sequence
    """
    def __init__(self, sequence: str, block: Block):
        # self.logger = Logger()
        self.sequence = sequence
        self.block = block
        self.__Clean()
        self.__Validate()
        return

    def __str__(self):
        return self.sequence

    def __Clean(self):
        self.sequence = self.sequence.strip().upper()
        return

    def __Validate(self):
        for block in self.sequence:
            if block not in self.block:
                raise InvalidSequenceError(f"InvalidSequenceError: {block} is not valid")
        return

class Score:
    """Scoring scheme for pairwise alignment

    Attributes:
        matrix          -- Scoring matrix
        existence       -- Value for gap existence
        extension       -- Value for gap extension
        block           -- Building block of sequence

    Methods:
        Matrix          -- Prints scoring matrix according to format indicated
        ValidateInput   -- Validates all input and automatically sets block type
    """
    def __init__(self, matrix: dict, block: Block, existence: int, extension: int):
        # self.logger = Logger()
        self.matrix = matrix
        self.existence = existence
        self.extension = extension
        self.block = block

        self.__Validate()

        return

    def __str__(self):
        """
        Output Format
        =====================
        | GapExistence: -11 |
        | GapExtension: -1  |
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
        times_across = len(self.block) + 1
        output = "="*(4*times_across+1) + f"\n| GapExistence: {str(self.existence).rjust(3).ljust(4)}|\n| GapExtension: {str(self.extension).rjust(3).ljust(4)}|\n"
        output += "="*(4*times_across+1) + "\n|   |"
        for hblock in self.block:
            output += f"{hblock.rjust(2).ljust(3)}|"
        for vblock in self.block:
            output += "\n" + "|---"*times_across + "|\n" + f"|{vblock.rjust(2).ljust(3)}|"
            for hblock in self.block:
                output += f"{str(self.matrix[vblock][hblock]).rjust(2).ljust(3)}|"
        output += "\n" + "="*(4*times_across+1)
        return output

    def __Validate(self):
        size = len(self.matrix.keys())
        if size != len(self.block):
            raise InvalidMatrixError(f"InvalidMatrixError: Size Given = {size} | Expected = {len(self.block)}")

        for key in self.matrix.keys():
            if key not in self.block:
                raise InvalidMatrixError(f"InvalidMatrixError: {key} is not a valid amino acid")

        for block in self.block:
            if block not in self.matrix.keys():
                raise InvalidMatrixError(f"InvalidMatrixError: {block} is not found in matrix")

        self.existence = -abs(self.existence)
        self.extension = -abs(self.extension)
        return

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
    def __init__(self, sequence1: str, sequence2: str, score: Score):
        # self.logger = Logger()
        self.score = score
        self.hseq = Sequence(sequence1, self.score.block)
        self.vseq = Sequence(sequence2, self.score.block)
        self.dp_array = None
        self.direction = None
        self.paths = None
        self.alignments = None

        self.__matrices = {
            "h_gap": np.zeros((len(self.vseq.sequence)+1, len(self.hseq.sequence)+1)),
            "v_gap": np.zeros((len(self.vseq.sequence)+1, len(self.hseq.sequence)+1)),
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
        self.dp_array = np.zeros((vlen, hlen))
        self.direction = [["" for i in range(hlen)] for j in range(vlen)]

        # Initialize recurrence matrices
        self.__matrices["v_gap"][0][0] = -np.Inf
        for i in range(1, hlen):
            self.__matrices["v_gap"][0][i] = exist + i*extend
            self.dp_array[0][i] = exist + i*extend
            self.direction[0][i] += "L"

        self.__matrices["h_gap"][0][0] = -np.Inf
        for j in range(1, vlen):
            self.__matrices["h_gap"][j][0] = exist + j*extend
            self.dp_array[j][0] = exist + j*extend
            self.direction[j][0] += "U"

        self.__matrices["match"][0][0] = -np.Inf

        # Affine gap penalty recurrence relation
        for j in range(1, vlen):
            for i in range(1, hlen):
                self.__matrices["h_gap"][j][i] = max(self.__matrices["h_gap"][j][i-1] + extend,
                                         self.dp_array[j][i-1] + exist + extend)
                self.__matrices["v_gap"][j][i] = max(self.__matrices["v_gap"][j-1][i] + extend,
                                         self.dp_array[j-1][i] + exist + extend)
                self.__matrices["match"][j][i] = (self.dp_array[j-1][i-1]
                                           + matrix[vseq[j-1]][hseq[i-1]])
                self.dp_array[j][i] = max(self.__matrices["h_gap"][j][i],
                                          self.__matrices["v_gap"][j][i],
                                          self.__matrices["match"][j][i])
                # Direction
                if self.dp_array[j][i] == self.__matrices["h_gap"][j][i]:
                    self.direction[j][i] += "L"
                if self.dp_array[j][i] == self.__matrices["v_gap"][j][i]:
                    self.direction[j][i] += "U"
                if self.dp_array[j][i] == self.__matrices["match"][j][i]:
                    self.direction[j][i] += "D"

        # Find alignment(s)
        self.direction = np.array(self.direction)
        self.alignments = []
        self.paths = []
        self.__GlobalPaths("", [], hlen-1, vlen-1)
        for i in range(len(self.alignments)):
            print(f"Alignment {i+1}\n{self.alignments[i]}")
        print(f"Score: {self.dp_array[vlen-1][hlen-1]}")
        return

    def __GlobalPaths(self, temp_path: str, nodes: list, i: int, j: int):
        """This function was rough to code...

        Args:
            temp_path (list): Temporary path
            nodes (list): Traversed nodes
            i (int): Optimal horizontal index in dp_array
            j (int): Optimal vertical index in dp_array
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
        self.dp_array = np.zeros((vlen, hlen))
        self.direction = [["" for i in range(hlen)] for j in range(vlen)]

        # Affine gap penalty recurrence relation
        for j in range(1, vlen):
            for i in range(1, hlen):
                self.__matrices["h_gap"][j][i] = max(self.__matrices["h_gap"][j][i-1] + extend,
                                         self.dp_array[j][i-1] + exist + extend, 0)
                self.__matrices["v_gap"][j][i] = max(self.__matrices["v_gap"][j-1][i] + extend,
                                         self.dp_array[j-1][i] + exist + extend, 0)
                self.__matrices["match"][j][i] = (self.dp_array[j-1][i-1]
                                           + matrix[vseq[j-1]][hseq[i-1]])
                self.dp_array[j][i] = max(self.__matrices["v_gap"][j][i],
                                          self.__matrices["h_gap"][j][i],
                                          self.__matrices["match"][j][i], 0)
                # Direction
                if self.dp_array[j][i] == self.__matrices["h_gap"][j][i]:
                    self.direction[j][i] += "L"
                if self.dp_array[j][i] == self.__matrices["v_gap"][j][i]:
                    self.direction[j][i] += "U"
                if self.dp_array[j][i] == self.__matrices["match"][j][i]:
                    self.direction[j][i] += "D"

        # Create alignment(s)
        self.direction = np.array(self.direction)
        self.alignments = []
        self.paths = []
        optimum = np.max(self.dp_array)
        for j in range(vlen):
            for i in range(hlen):
                if self.dp_array[j][i] == optimum:
                    self.__LocalPaths("", [], i, j)
        for i in range(len(self.alignments)):
            print(f"Alignment {i+1}\n{self.alignments[i]}")
        print(f"Score: {optimum}")
        return

    def __LocalPaths(self, temp_path: str, nodes: list, i: int, j: int):
        """This one was easier...

        Args:
            temp_path (list): Temporary path
            i (int): Optimal horizontal index in dp_array
            j (int): Optimal vertical index in dp_array
        """
        # Perform depth first walk to nearest 0 while obtaining matching
        if self.dp_array[j][i] != 0:
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
        axis[0].table(cellText = self.dp_array,
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
