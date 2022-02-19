from Error import InvalidSequenceError, InvalidSequenceTypeError

class Sequence:
    """_summary_
    """
    def __init__(self, sequence: str, sequenceType):
        self.sequence = sequence
        self.sequenceType = sequenceType
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
        if self.sequenceType not in [NUCLEOTIDES, AMINO_ACIDS]:
            raise InvalidSequenceTypeError()
        
        for monomer in self.sequence:
            if monomer not in self.sequenceType:
                raise InvalidSequenceError(f"InvalidSequenceError: {monomer} is not valid")
        return
    
    def __Transform(self):
        if self.sequenceType == AMINO_ACIDS:
            self.__class__ = AASequence
        if self.sequenceType == NUCLEOTIDES:
            self.__class__ = NTSequence
        return
            
    def FindSubsequence(self, subsequence: str, overlapping=True):
        """Find all occurrences of subsequences in forward strand.

        Args:
            subsequence (_type_): _description_
            overlapping (bool, optional): _description_. Defaults to True.

        Raises:
            InvalidSequenceError: _description_
        """
        for monomer in subsequence:
            if monomer not in self.sequenceType:
                raise InvalidSequenceError(f"InvalidSequenceError: {monomer} is not valid")
    
        indices = []
        index = self.sequence.find(subsequence)
        while index != -1:
            indices.append(index)
            index = self.sequence.find(subsequence,
                                       index + (1 if overlapping else len(subsequence)))
        return
        
class AASequence(Sequence):
    
    def ToDNASequence(self):
        return
        
class NTSequence(Sequence):
    
    def Complement(self):
        """Computes the reverse complement of a sequence.

        Returns:
            string: Reverse complement of a sequence from 5' to 3' direction.
        """
        complement = ""
        for monomer in self.sequence:
            complement += COMPLEMENT[monomer]
            
        return NTSequence(complement[::-1])
    
    def FindFRSubsequence(self, subsequence: str, overlapping=True):
        """Finds all occurrences of subsequence in forward and reverse strand

        Args:
            subsequence (str): _description_
            overlapping (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """
        forwardIndices = self.FindSubsequence(subsequence, overlapping)
        reverseIndices = self.Complement().FindSubsequence(subsequence, overlapping)
        indices = {"forward": forwardIndices, "reverse": reverseIndices}
        return indices
    
    def DNA2RNA(self):
        sequence = self.sequence
        sequence.replace('T', 'U')
        return sequence
    
    def ToAASequence(self, sequence):
        translatedSequence = ""
        for i in range(0, len(sequence), 3):
            translatedSequence += CODON[sequence[i:i+3]]
        return AASequence(translatedSequence)
    
    def To6AASequences(self):
        """_summary_
        """
        AASequences = {"forward": [], "reverse": []}
        reverseSequence = self.Complement().sequence
        for i in range(0, 3):
            AASequences["forward"].append(self.ToAASequence(self.sequence[i:]))
            AASequences["reverse"].append(self.ToAASequence(reverseSequence[i:]))
        return AASequences
    
    

# Constants
NUCLEOTIDES = ['A', 'C', 'G', 'T']
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
COMPLEMENT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

CODON = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}
