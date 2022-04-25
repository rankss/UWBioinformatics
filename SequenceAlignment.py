import numpy as np
from dataclasses import dataclass
from typing import Literal
from Sequence import Sequence
from Score import Score
from Cluster import Cluster

class PWA:
    """Represents a Global/Local Pairwise Alignment
    """
    @dataclass
    class PWAData:
        score: float
        alignments: list[tuple]
        
        def distance(self):
            minDistance = np.inf
            for alignment in self.alignments:
                score = 0.0
                for top, bot in zip(alignment[0], alignment[1]):
                    if top == '-' or bot == '-':
                        score += 0.75
                    elif top != bot:
                        score += 1
                minDistance = min(score, minDistance)
            
            return minDistance
                
    _GLOBAL = 0
    _LOCAL = 1
    _DIR_DICT = {
        "L": [-1, 0], # Left
        "U": [0, -1], # Up
        "D": [-1, -1] # Diagonal
    }
    
    def _path(hSeq: Sequence, vSeq: Sequence, temp_path: str, dpArray: np.ndarray, dirArray: list[list[str]], col: int, row: int, data: PWAData, alignType: Literal=_GLOBAL) -> None:
        if alignType == PWA._GLOBAL:
            comparison = col > 0 or row > 0
        if alignType == PWA._LOCAL:
            comparison = dpArray[row, col] != 0
            
        alignments = []
            
        if comparison:
            curr = dirArray[row][col]
            for direction in curr:
                PWA._path(hSeq, vSeq, temp_path + direction,
                           dpArray, dirArray,
                           col + PWA._DIR_DICT[direction][0],
                           row + PWA._DIR_DICT[direction][1],
                           data, alignType)
        else:
            temp_path = temp_path[::-1]
            hMatch, vMatch = "", ""
            hCounter, vCounter = 0, 0
            for direction in temp_path:
                if direction == "D":
                    hMatch += hSeq[hCounter]
                    vMatch += vSeq[vCounter]
                    hCounter += 1
                    vCounter += 1
                if direction == "L":
                    hMatch += hSeq[hCounter]
                    vMatch += "-"
                    hCounter += 1
                if direction == "U":
                    hMatch += "-"
                    vMatch += vSeq[vCounter]
                    vCounter += 1
            match = (hMatch, vMatch)
            data.alignments.append(match)
        return
        
    @staticmethod
    def Global(hSeq: Sequence, vSeq: Sequence, score: Score) -> PWAData:
        hLen, vLen = len(hSeq) + 1, len(vSeq) + 1
        exist, extend = score.existence, score.extension
        matrix = score.matrix
        
        # Define recurrence matrices
        dpArray = np.zeros((vLen, hLen))
        direction = [["" for i in range(hLen)] for j in range(vLen)]
        vGap = np.full((vLen, hLen + 1), -np.inf)
        hGap = np.full((vLen, hLen + 1), -np.inf)
        match = np.full((vLen, hLen + 1), -np.inf)
        
        for col in range(1, hLen):
            vGap[0, col] = exist + col*extend
            dpArray[0, col] = exist + col*extend
            direction[0][col] += "L"

        for row in range(1, vLen):
            hGap[row, 0] = exist + row*extend
            dpArray[row, 0] = exist + row*extend
            direction[row][0] += "U"

        # Affine gap penalty recurrence relation
        for row in range(1, vLen):
            for col in range(1, hLen):
                hGap[row, col] = max(hGap[row, col-1] + extend,
                                     dpArray[row, col-1] + exist + extend)
                
                vGap[row, col] = max(vGap[row-1, col] + extend,
                                     dpArray[row-1, col] + exist + extend)
                
                match[row, col] = dpArray[row-1, col-1] + matrix[vSeq[row-1]][hSeq[col-1]]
                
                dpArray[row, col] = max(hGap[row, col], vGap[row, col], match[row, col])
                
                if dpArray[row, col] == hGap[row, col]:
                    direction[row][col] += "L"
                if dpArray[row, col] == vGap[row, col]:
                    direction[row][col] += "U"
                if dpArray[row, col] == match[row, col]:
                    direction[row][col] += "D"

        optimal = dpArray[vLen-1][hLen-1]
        data = PWA.PWAData(optimal, [])
        PWA._path(hSeq, vSeq, "", dpArray, direction, hLen-1, vLen-1, data, PWA._GLOBAL)
        return data

    @staticmethod
    def Local(hSeq: Sequence, vSeq: Sequence, score: Score) -> PWAData:
        hLen, vLen = len(hSeq) + 1, len(vSeq) + 1
        exist, extend = score.existence, score.extension
        matrix = score.matrix
        
        # Define recurrence matrices
        dpArray = np.zeros((vLen, hLen))
        vGap = np.zeros((vLen, hLen + 1))
        hGap = np.zeros((vLen, hLen + 1))
        match = np.zeros((vLen, hLen + 1))
        
        # Affine gap penalty recurrence relation
        for row in range(1, vLen):
            for col in range(1, hLen):
                hGap[row, col] = max(0, hGap[row, col-1] + extend,
                                     dpArray[row, col-1] + exist + extend)
                
                vGap[row, col] = max(0, vGap[row-1, col] + extend,
                                     dpArray[row-1, col] + exist + extend)
                
                match[row, col] = dpArray[row-1, col-1] + matrix[vSeq[row-1]][hSeq[col-1]]
                
                dpArray[row, col] = max(0, hGap[row, col], vGap[row, col], match[row, col])
                
        # Create alignment(s)
        alignments = []
        paths = []
        optimal = np.max(dpArray)
        for row in range(vLen):
            for col in range(hLen):
                if dpArray[row, col] == optimal:
                    PWA._path("", col, row, PWA._LOCAL)
            
        return

class MSA:
    """Represents a multiple sequence alignment.
    """
    @dataclass
    class MSAData:
        pass
    
    @staticmethod
    def clustalw(sequences: list[Sequence], score: Score):
        # Produce distance matrix
        distMatrix = np.zeros((len(sequences), len(sequences)))
        for row, vSeq in enumerate(sequences[:-1]):
            for col, hSeq in sequences[row+1:]:
                distance = PWA.Global(hSeq, vSeq, score).distance()
                distMatrix[row, col] = distance
                distMatrix[col, row] = distance
        print(distMatrix)

        # Produce guide tree
        cluster = Cluster(distMatrix, [sequence.taxa for sequence in sequences])
        digraph = cluster.upgma()
        
        # Follow guide tree for alignment
        return