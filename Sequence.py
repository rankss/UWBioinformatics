from Error import InvalidSequenceError, InvalidSequenceTypeError
from typing import Literal
from dataclasses import dataclass

@dataclass
class Sequence:
    """_summary_
    """
    
    NUCLEOTIDES = ['A', 'C', 'G', 'T']
    AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    COMPLEMENT = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    CODON_DICT = {
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
    
    def __init__(self, sequence: str, taxa: str, sequenceType: Literal=None):
        self.sequence = sequence
        self.taxa = taxa
        self.sequenceType = sequenceType
        self.summary = {}
        
        self._clean()
        if self.sequenceType is None:
            self._autoDetectSequenceType()
        self._validate()
        if type(self) is Sequence:
            self._transform()
        return

    def __str__(self):
        return self.sequence
    
    def __getitem__(self, i):
        return self.sequence[i]
    
    def __len__(self):
        return len(self.sequence)
    
    def _clean(self):
        self.sequence = self.sequence.strip().upper()
        return
    
    def _validate(self):
        if self.sequenceType not in [Sequence.NUCLEOTIDES, Sequence.AMINO_ACIDS]:
            raise InvalidSequenceTypeError("Sequence type is not valid.")
    
    def _transform(self):
        if self.sequenceType == Sequence.NUCLEOTIDES:
            self.__class__ = NTSequence
            return
        if self.sequenceType == Sequence.AMINO_ACIDS:
            self.__class__ = AASequence
            return
        return
    
    def _autoDetectSequenceType(self):
        monomers = set(self.sequence)
        if monomers.issubset(set(Sequence.NUCLEOTIDES)):
            self.sequenceType = Sequence.NUCLEOTIDES
            return
        if monomers.issubset(set(Sequence.AMINO_ACIDS)):
            self.sequenceType = Sequence.AMINO_ACIDS
            return
        raise InvalidSequenceError("Sequence is neither a protein nor DNA.")
            
    def findSubsequence(self, subsequence: str, overlapping=True):
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
        return indices
    
    def summary(self):
        self.summary["length"] = len(self.sequence)
        
        self.summary["frequency"] = {monomer:0 for monomer in self.sequenceType}
        for monomer in self.sequence:
            self.summary["frequency"][monomer] += 1
        return
        
class AASequence(Sequence):
    
    def __init__(self, sequence: str, taxa: str, sequenceType: Literal=Sequence.AMINO_ACIDS):
        super().__init__(sequence, taxa, sequenceType)
    
    def toDNASequence(self):
        return
        
class NTSequence(Sequence):
    
    def __init__(self, sequence: str, taxa: str, sequenceType=Sequence.NUCLEOTIDES):
        super().__init__(sequence, taxa, sequenceType)
    
    def complement(self):
        """Computes the reverse complement of a sequence.

        Returns:
            string: Reverse complement of a sequence from 5' to 3' direction.
        """
        complement = ""
        for monomer in self.sequence:
            complement += Sequence.COMPLEMENT[monomer]
            
        return NTSequence(complement[::-1], f"{self.taxa}_complement")
    
    def findFRSubsequence(self, subsequence: str, overlapping=True):
        """Finds all occurrences of subsequence in forward and reverse strand

        Args:
            subsequence (str): _description_
            overlapping (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """
        indices = {
            "forward": self.findSubsequence(subsequence, overlapping), 
            "reverse": self.complement().findSubsequence(subsequence, overlapping)
        }
        
        return indices
    
    def to3AASequences(self):
        AASequences = []
        for i in range(3):
            translatedSequence = ""
            for j in range(i, len(self.sequence), 3):
                codon = self.sequence[j:j+3]
                if codon not in Sequence.CODON_DICT.keys():
                    break
                if Sequence.CODON_DICT[codon] == "_":
                    break
                translatedSequence += Sequence.CODON_DICT[codon]
            AASequences.append(AASequence(translatedSequence, f"{self.taxa}_translated{i+1}"))
        return AASequences
    
    def toFR3AASequences(self):
        """_summary_
        """
        reverseSequence = self.complement()
        AASequences = {
            "forward": self.toAASequence(),
            "reverse": reverseSequence.toAASequence()
        }
        
        return AASequences
    
    def tRNAScan(self):
        """Performs tRNA decision tree
        """
        windowSize = 76
        score = 0
        for i in range(0, len(self.sequence) - windowSize + 1):
            window = self.sequence[i:i+windowSize]
            # Count T-phi-C signal invariant bases
            # count = 0
            # invariantBases = {52: 'G', 54: 'T', 55: 'C', 60: 'C'}
            # for key in invariantBases.keys():
            #     if window[key] == invariantBases[key]:
            #         count += 1
            # if count > 2:
            #     score += 1
            # else:
            #     break
            # # Check T-phi-C arm base pairing
            # count = 0
            # basepairing = {52: 60, 51: 61, 50: 62, 49: 63, 48: 64}
            # for key in basepairing.keys():
            #     if COMPLEMENT[window[key]] == window[basepairing[key]]:
            #         count += 1
            # if count > 4:
            #     score += 1
            # else:
            #     break
            # # Count D signal invariant bases
            # count = 0
            # invariantBases = {7: 'T', 9: 'G', 13: 'A'}
            # for key in invariantBases.keys():
            #     if window[key] == invariantBases[key]:
            #         count += 1
            # if count > 2:
            #     score += 1
            # else:
            #     break
            # # Check D arm base pairing
            # count = 0
            # basepairing = {9: 24, 10: 23, 11: 22}
            # for key in basepairing.keys():
            #     if COMPLEMENT[window[key]] == window[basepairing[key]]:
            #         count += 1
            # if count > 2:
            #     score += 1
            # else:
            #     break
        return
    
    def TM(self):
        pass
    
# Constants
    



