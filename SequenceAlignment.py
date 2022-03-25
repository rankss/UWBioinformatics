import numpy as np
from dataclasses import dataclass
from typing import Literal
from Sequence import Sequence
from Score import Score
from Cluster import Cluster

class PWA:
    """Global/Local Pairwise Alignment
    """
    @dataclass
    class PWAData:
        score: float
        paths: list
        alignments: list
    
    GLOBAL = 0
    LOCAL = 1
    __DIR_DICT = {
        "L": [-1, 0], # Left
        "U": [0, -1], # Up
        "D": [-1, -1] # Diagonal
    }

    # @staticmethod
    # def __match(self, temp_path: str): 
    #     # Charlotte TODO Match has bug for local alignment
    #     hMatch, vMatch = "", ""
    #     hCounter, vCounter = 0, 0
    #     for direction in temp_path:
    #         if direction == "D":
    #             hMatch += self.hSeq[hCounter]
    #             vMatch += self.vSeq[vCounter]
    #             hCounter += 1
    #             vCounter += 1
    #         if direction == "L":
    #             hMatch += self.hSeq[hCounter]
    #             vMatch += "-"
    #             hCounter += 1
    #         if direction == "U":
    #             hMatch += "-"
    #             vMatch += self.vSeq[vCounter]
    #             vCounter += 1
    #     match = f"{hMatch}\n{vMatch}"
    #     self.alignments.append(match)
    #     return
    
    def __path(temp_path: str, nodes: list, dpArray: np.ndarray, direction: list, col: int, row: int, alignment: Literal=GLOBAL) -> tuple:
        if alignment == PWA.GLOBAL:
            comparison = col > 0 or row > 0
        if alignment == PWA.LOCAL:
            comparison = dpArray[row, col] != 0
            
        if comparison:
            curr = direction[row][col]
            for direction in curr:
                PWA.__path(temp_path + direction, nodes + [(col, row)],
                           dpArray, direction,
                           col + PWA.__DIR_DICT[direction][0],
                           row + PWA.__DIR_DICT[direction][1],
                           alignment)
        
        # TODO Charlotte: Insert match code here
        # 
        
        return nodes, PWA.__match(temp_path[::-1])
    
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
        PWA.__path("", [], hLen-1, vLen-1, PWA.GLOBAL)
        
        return

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
                    PWA.__path("", [], col, row, PWA.LOCAL)
                
        return

class MSA:
    """[summary]
    """
    @dataclass
    class MSAData:
        pass
    
    def clustalw(sequences: list, score: Score):
        # Produce distance matrix
        distMatrix = np.zeros((len(sequences), len(sequences)))
        for row, vSeq in enumerate(sequences[:-1]):
            for col, hSeq in enumerate(sequences[row+1:], row+1):
                distance = PWA.Global(hSeq, vSeq, score)

                
        # Produce guide tree
        cluster = Cluster(distMatrix, [sequence.taxa for sequence in sequences])
        root = cluster.upgma()
        
        # Follow guide tree for alignment
        return