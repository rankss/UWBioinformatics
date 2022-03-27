import numpy as np
from itertools import product
from Graph import Node, Digraph

class Cluster:
    def __init__(self, distMatrix: np.ndarray, taxa: list):
        self.distMatrix = np.array(distMatrix, dtype=float)
        np.fill_diagonal(self.distMatrix, np.inf)
        self.taxa = taxa
        
    def upgma(self) -> Digraph:
        """Performs UPGMA on a distance matrix with corresponding taxa.
        """
        distMatrix = self.distMatrix
        # Initialize list of node and maps of nodes -> index
        digraph = Digraph([Node(taxon) for taxon in self.taxa])
        nodeToIndex = {node:i for i, node in enumerate(digraph.tmpNodes)}
        
        for i in range(1, len(self.taxa)):
            # Find min of distMatrix and obtain their corresponding nodes and assign distance
            index = np.unravel_index(distMatrix.argmin(), distMatrix.shape)
            distance = distMatrix[index[0], index[1]]
            rowNode, colNode = digraph.tmpNodes[index[0]], digraph.tmpNodes[index[1]]
            
            # Merge row/col Nodes into new node, remove them and append merged node
            rowDist = distance/2 - (digraph.adjList[rowNode][0].distance if len(digraph.adjList[rowNode]) else 0)
            colDist = distance/2 - (digraph.adjList[colNode][0].distance if len(digraph.adjList[colNode]) else 0)
            mergeNode = digraph.join((rowNode, colNode), (rowDist, colDist))
            
            # Reinitialize list of node and maps of nodes -> index
            newMatrix = np.full((len(self.taxa)-i, len(self.taxa)-i), np.inf)
            newNodeToIndex = {node:i for i, node in enumerate(digraph.tmpNodes)}
                
            # Loop in upper triangle fashion, lower triangle can be filled by swapping indices
            for j, rowNode in enumerate(digraph.tmpNodes[:-1]):
                for colNode in digraph.tmpNodes[j+1:]:
                    if colNode != mergeNode:
                        newMatrix[newNodeToIndex[rowNode], newNodeToIndex[colNode]] = distMatrix[nodeToIndex[rowNode], nodeToIndex[colNode]]
                        newMatrix[newNodeToIndex[colNode], newNodeToIndex[rowNode]] = distMatrix[nodeToIndex[colNode], nodeToIndex[rowNode]]
                    else:
                        nodeProduct = product([rowNode], digraph.adjList[colNode])
                        totalDistance = sum([distMatrix[nodeToIndex[node1], nodeToIndex[node2]]*node2.leaves for node1, node2 in nodeProduct])
                        newMatrix[newNodeToIndex[rowNode], newNodeToIndex[colNode]] = totalDistance/colNode.leaves
                        newMatrix[newNodeToIndex[colNode], newNodeToIndex[rowNode]] = totalDistance/colNode.leaves

            # Update all
            distMatrix = newMatrix
            nodeToIndex = newNodeToIndex
        
        return digraph
    
    def neighborJoining(self) -> Node:
        # Charlotte TODO
        
        pass


    