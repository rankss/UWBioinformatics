from Constants import Block, Translation
from Error import InvalidSequenceError

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
        if type(self) == Sequence:
            self.__Transform()
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
    
    def __Transform(self):
        if self.block == Block.AMINO_ACIDS:
            self.__class__ = AASequence
        else:
            self.__class__ = NTSequence
            
    def FindSubsequence(self, subsequence):
        for block in subsequence:
            if block not in self.block:
                raise InvalidSequenceError(f"InvalidSequenceError: {block} is not valid")
    
        index = []
        
class AASequence(Sequence):
    
    pass
        
class NTSequence(Sequence):
    
    def Complement(self):
        """Computes the complement of nucleotide sequence
        """
        complement = ""
        for block in self.sequence:
            complement += COMPLEMENT[block]
            
        return complement[::-1]
    
    def ToAASequence(self):
        """Computes the 6 ORF of nucleotide sequence
        """
        