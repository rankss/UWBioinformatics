from Constants import *
from Logger import Logger
from Error import *

class Sequence:
    """
    Base class for sequence strings

    Attributes:
        sequence    -- Stores sequence information after stripping and converting to upper case

    Methods:
        Validate    -- Implemented by subclasses
        Clean       -- Basic cleaning of sequence
    """
    def __init__(self, sequence: str):
        self.logger = Logger()
        self.sequence = sequence
        self.Clean()
        return

    def __str__(self):
        return self.sequence

    def Clean(self):
        self.sequence = self.sequence.strip().upper()
        return

    def Validate(self):
        raise NotImplementedError("Implemented by subclasses")

class DNASequence(Sequence):
    """
    Methods:
        Validate -- Checks if nucleotide sequence string is valid
    """
    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.Validate()

    def Validate(self):
        for nucleotide in self.sequence:
            if nucleotide not in BlockConstants.NUCLEOTIDES:
                raise InvalidSequenceError(f"InvalidSequenceError: {nucleotide} is not a nucleotide")
        return

class AASequence(Sequence):
    """
    Methods:
        Validate -- Checks if amino acid sequence string is valid
    """
    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.Validate()

    def Validate(self):
        for amino_acid in self.sequence:
            if amino_acid not in BlockConstants.AMINO_ACIDS:
                raise InvalidSequenceError(f"InvalidSequenceError: {amino_acid} is not an amino acid")
        return
