from Sequence import Sequence
from Graph import Digraph

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
    
    @staticmethod
    def Newick(newick: str) -> Digraph:
        """Converts newick string into rooted tree.
        """
        def findComma(newick: str) -> list:
            """Finds all comma indices between first and last bracket appended with -1.
            """
            bracketCount = 0
            commaList = []
            for i, ch in enumerate(newick):
                if ch == '(':
                    bracketCount += 1
                if ch == ')':
                    bracketCount -= 1
                if ch == ',' and bracketCount == 1:
                    commaList.append(i)
            commaList.append(-1)
            return commaList
    