from Sequence import Sequence

class Parser:
    
    @staticmethod
    def Fasta(filename: str):
        with open(filename) as file:
            lines = [line.strip() for line in file]
        
        sequence = ""
        sequenceCollection = []
        for line in lines:
            if line[0] == '>':
                if len(sequence):
                    sequenceCollection.append(Sequence(sequence, sequenceName))
                    sequence = ""
                sequenceName = line[1:]
            else:
                sequence += line
        sequenceCollection.append(Sequence(sequence, sequenceName))
        
        return sequenceCollection
    