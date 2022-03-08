import numpy as np
from itertools import product

class Node:
    def __init__(self, taxa: str=None, distance: float=0.0):
        self.taxa = taxa
        self.distance = distance
        if type(self.taxa) is str:
            self.isLeaf = True
        else:
            self.isLeaf = False
            
    def __len__(self) -> int:
        if self.isLeaf:
            return 1
        length = 0
        for leaf in self.taxa:
            length += len(leaf)
        return length
    
    def updateDistance(self, distance: float):
        if self.isLeaf:
            self.distance = distance
        else:
            self.distance = distance - max([node.totalDistance() for node in self.taxa])
        return
    
    def children(self) -> list:
        if self.isLeaf:
            return [self]
        return self.taxa
    
    def totalDistance(self) -> float:
        if self.isLeaf:
            return self.distance
        return self.distance + max([node.totalDistance() for node in self.taxa])
    
    def toNewick(self):
        if self.isLeaf:
            return f"{self.taxa}:{self.distance}"
        else:
            return f"({','.join([node.toNewick() for node in self.taxa])}):{self.distance}"
    
    @staticmethod
    def newickToNode(newick: str):
        pass

class Cluster:
    def __init__(self, distMatrix: np.ndarray, taxa: list):
        self.distMatrix = np.array(distMatrix, dtype=float)
        np.fill_diagonal(self.distMatrix, np.inf)
        self.taxa = taxa
        
    def upgma(self) -> Node: 
        # https://codereview.stackexchange.com/questions/263416/upgma-tree-building-in-python with a few modifications
        distMatrix = self.distMatrix
        nodes = [Node(taxon) for taxon in self.taxa] # list of initialized nodes
        nodeToIndex = {node:i for i, node in enumerate(nodes)} # define maps node -> index
        indexToNode = {i:node for i, node in enumerate(nodes)} # define maps index -> node
        
        iterations = len(self.taxa)
        for i in range(1, iterations): # UPGMA always run for n - 1 iterations, where n = length of taxa
            index = np.unravel_index(distMatrix.argmin(), distMatrix.shape) # finds min index of ndarray
            value = distMatrix[index[0], index[1]] # obtains minimum value
            rowNode, colNode = indexToNode[index[0]], indexToNode[index[1]] # obtain nodes from indices
            rowNode.updateDistance(value/2), colNode.updateDistance(value/2)
                        
            mergeNode = Node([rowNode, colNode]) # merge rowNode and colNode
            nodes.remove(rowNode), nodes.remove(colNode) # remove rowNode and colNode
            nodes.append(mergeNode) # append mergeNode to list of nodes
            
            newMatrix = np.full((iterations-i, iterations-i), np.inf) # define new distMatrix of infinity
            newNodeToIndex = {node:i for i, node in enumerate(nodes)} # define new map node -> index
            newIndexToNode = {i:node for i, node in enumerate(nodes)} # define new map index -> node
            
            # Loop in upper triangle fashion, lower triangle can be filled by swapping indices
            for i, rowNode in enumerate(nodes[:-1]):
                newRowIndex = newNodeToIndex[rowNode]
                rowChild = [rowNode] # we actually only want the rowNode itself here
                for colNode in nodes[i+1:]:
                    newColIndex = newNodeToIndex[colNode]
                    colChild = colNode.children() if colNode == mergeNode else [colNode] # if node is mergeNode, we want its children for arithmetic mean
                    nodeProduct = list(product(rowChild, colChild)) # cartesian product of two lists

                    # Black magic calculation
                    result = np.array([[distMatrix[nodeToIndex[rowNode], nodeToIndex[colNode]]*len(colNode), len(colNode)] for rowNode, colNode in nodeProduct])
                    newMatrix[newRowIndex, newColIndex] = result[:,0].sum()/result[:,1].sum()
                    newMatrix[newColIndex, newRowIndex] = result[:,0].sum()/result[:,1].sum()
            
            # Update all
            distMatrix = newMatrix
            nodeToIndex = newNodeToIndex
            indexToNode = newIndexToNode
        
        return nodes[0]
    
    def neighborJoining(self) -> Node:
        pass

class Newick:
    @staticmethod
    def toNewick(root: Node):
        if root.isLeaf:
            return f"{root.taxa}:{root.distance}"
        else:
            return f"({','.join([Newick.toNewick(node) for node in root.taxa])}):{root.distance}"
    
    @staticmethod
    def toTree(newick: str):
        pass
        