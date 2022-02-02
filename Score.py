from Constants import BlockConstants
from Logger import Logger
from Error import *

class Score:
    """
    Scoring scheme for pairwise alignment

    Attributes:
        matrix          -- Scoring matrix
        existence       -- Value for gap existence
        extension       -- Value for gap extension
        block           -- Building block of sequence

    Methods:
        Matrix          -- Prints scoring matrix according to format indicated
        ValidateInput   -- Validates all input and automatically sets block type
    """
    def __init__(self, matrix: dict, block: BlockConstants, existence: int, extension: int):
        # self.logger = Logger()
        self.matrix = matrix
        self.existence = existence
        self.extension = extension
        self.block = block

        self.Validate()

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

    def Validate(self):
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
        