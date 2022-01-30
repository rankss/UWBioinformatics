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
        nucleotides = ["A", "G", "C", "T"]
        for nt in self.sequence:
            if nt not in nucleotides:
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
        aminoacids = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", 
                        "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "T"]
        for aa in self.sequence:
            if aa not in aminoacids:
                raise SequenceInvalidError(self.sequence)
        return
