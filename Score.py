import Constants
from Logger import Logger

class Score:
    """
    Scoring scheme for pairwise alignment

    Attributes:
        matrix          -- Scoring matrix
        existence       -- Value for gap existence (maximum = 99)
        extension       -- Value for gap extension (maximum = 99)
        block           -- Building block of sequence

    Methods:
        Matrix          -- Prints scoring matrix according to format indicated
        ValidateInput   -- Validates all input and automatically sets block type
    """
    def __init__(self, matrix: dict, existence: int, extension: int):
        self.matrix = matrix
        self.existence = existence
        self.extension = extension
        self.block = None

        self.ValidateInput()

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
        timesAcross = len(sequence) + 1
        output = f"="*(4*timesAcross+1) + f"\n| GapExistence: {str(self.existence).rjust(3).ljust(4)}|\n| GapExtension: {str(self.extension).rjust(3).ljust(4)}|\n"
        output += f"="*(4*timesAcross+1) + f"\n|   |"
        for hnt in self.sequence:
            output += f"{hnt.rjust(2).ljust(3)}|"
        for vnt in self.sequence:
            output += f"\n" + f"|---"*timesAcross + f"|\n" + f"|{vnt.rjust(2).ljust(3)}|"
            for hnt in sequence:
                output += f"{str(self.matrix[vnt][hnt]).rjust(2).ljust(3)}|"
        output += f"\n" + f"="*((4*timesAcross+1)+1) + f"\n"
        return output

    def ValidateInput(self):
        # TODO
        return
        