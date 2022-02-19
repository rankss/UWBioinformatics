from Error import InvalidMatrixError, InvalidSequenceTypeError
from Sequence import NUCLEOTIDES, AMINO_ACIDS

class Score:
    """_summary_
    """
    def __init__(self, match: int, mismatch: int, existence: int, extension: int, sequenceType):
        self.matrix = None
        self.match = match
        self.mismatch = mismatch
        self.existence = existence
        self.extension = extension
        self.sequenceType = sequenceType
        self.__Construct()
        return

    def __str__(self):
        """Output Format
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
        times_across = len(self.sequenceType) + 1
        output = "="*(4*times_across+1) + f"\n| GapExistence: {str(self.existence).rjust(3).ljust(4)}|\n| GapExtension: {str(self.extension).rjust(3).ljust(4)}|\n"
        output += "="*(4*times_across+1) + "\n|   |"
        for hmonomer in self.sequenceType:
            output += f"{hmonomer.rjust(2).ljust(3)}|"
        for vmonomer in self.sequenceType:
            output += "\n" + "|---"*times_across + "|\n" + f"|{vmonomer.rjust(2).ljust(3)}|"
            for hmonomer in self.sequenceType:
                output += f"{str(self.matrix[vmonomer][hmonomer]).rjust(2).ljust(3)}|"
        output += "\n" + "="*(4*times_across+1)
        return output

    def __ValidateMatrix(self):
        size = len(self.matrix.keys())
        if size != len(self.sequenceType):
            raise InvalidMatrixError(f"InvalidMatrixError: Size Given = {size} | Expected = {len(self.sequenceType)}")

        for key in self.matrix.keys():
            if key not in self.sequenceType:
                raise InvalidMatrixError(f"InvalidMatrixError: {key} is not a valid amino acid")

        for monomer in self.sequenceType:
            if monomer not in self.matrix.keys():
                raise InvalidMatrixError(f"InvalidMatrixError: {monomer} is not found in matrix")
        return
    
    def __Construct(self):
        if self.sequenceType not in [NUCLEOTIDES, AMINO_ACIDS]:
            raise InvalidSequenceTypeError()

        self.matrix = {}
        for vmonomer in self.sequenceType:
            self.matrix[vmonomer] = {}
            for hmonomer in self.sequenceType:
                if vmonomer == hmonomer:
                    self.matrix[vmonomer][hmonomer] = self.match
                else:
                    self.matrix[vmonomer][hmonomer] = self.mismatch
                    
        self.existence = -abs(self.existence)
        self.extension = -abs(self.extension)
        return
    
    def Matrix(self, matrix: dict):
        self.matrix = matrix
        self.__ValidateMatrix()
        return
    
EXAMPLE_MATRIX = {
    'A': {'A': 1, 'G': 0, 'C': 0, 'T': 0},
    'G': {'A': 0, 'G': 1, 'C': 0, 'T': 0},
    'C': {'A': 0, 'G': 0, 'C': 1, 'T': 0},
    'T': {'A': 0, 'G': 0, 'C': 0, 'T': 1}
}

BLOSUM62 = {
    # TODO
}
