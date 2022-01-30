import Constants
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
        for nt in self.sequence:
            if nt not in Constants.NUCLEOTIDES:
                raise SequenceInvalidError(self.sequence)
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
        for aa in self.sequence:
            if aa not in Constants.AMINO_ACIDS:
                raise SequenceInvalidError(self.sequence)
        return
