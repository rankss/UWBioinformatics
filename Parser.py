from Sequence import Sequence

class Parser:
    
    @staticmethod
    def Fasta(filename: str):
        with open(filename, 'r') as file:
            lines = [line.strip() for line in file]
        
        sequence = ""
        taxa = ""
        sequenceCollection = []
        for line in lines:
            if line[0] == '>':
                if len(sequence):
                    sequenceCollection.append(Sequence(sequence, taxa))
                    sequence = ""
                taxa = line[1:]
            else:
                sequence += line
        sequenceCollection.append(Sequence(sequence, taxa))
        
        return sequenceCollection
    