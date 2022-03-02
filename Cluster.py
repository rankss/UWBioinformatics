import numpy as np

class Node:
    def __init__(self, taxa: list=[], distance: float=0.0):
        self.taxa = taxa
        self.distance = distance

class Tree:
    @staticmethod
    def Dendrogram(self, root: Node):
        pass

class Cluster:
    def __init__(self, distMatrix: np.ndarray, taxa: list):
        self.distMatrix = distMatrix
        self.taxa = taxa
        