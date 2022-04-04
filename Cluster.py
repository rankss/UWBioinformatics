import numpy as np
from itertools import product
from Graph import Node, Digraph

class Cluster:
    def __init__(self, distMatrix: np.ndarray, taxa: list[str]):
        self.distMatrix = np.array(distMatrix, dtype=float)
        self.taxa = taxa
        
    def upgma(self) -> Digraph:
        """Performs UPGMA on a distance matrix with corresponding taxa.
        """
        distMatrix = self.distMatrix
        np.fill_diagonal(distMatrix, np.inf)
        # Initialize list of node and maps of nodes -> index
        tmpNodes = [Node(taxon) for taxon in self.taxa]
        digraph = Digraph(tmpNodes)
        nodeToIndex = {node:i for i, node in enumerate(tmpNodes)}
        
        for i in range(1, len(self.taxa)):
            # Find min of distMatrix and obtain their corresponding nodes and assign distance
            row, col = np.unravel_index(distMatrix.argmin(), distMatrix.shape)
            distance = distMatrix[row, col]/2
            rowNode, colNode = tmpNodes[row], tmpNodes[col]
            rowNode.total, colNode.total = distance, distance
            rowNode.distance = distance - (digraph.adjList[rowNode][0].total if digraph.adjList[rowNode] else 0)
            colNode.distance = distance - (digraph.adjList[colNode][0].total if digraph.adjList[colNode] else 0)
            
            # Merge row/col Nodes into new node, remove them and append merged node
            mergeNode = digraph.join((rowNode, colNode), tmpNodes)
            print(f"Distance = {distance} | Merged: {rowNode} with {colNode}") # For debugging
            
            # Reinitialize list of node and maps of nodes -> index
            newMatrix = np.full((len(self.taxa)-i, len(self.taxa)-i), np.inf)
            newNodeToIndex = {node:i for i, node in enumerate(tmpNodes)}

            # Loop in upper triangle fashion, lower triangle can be filled by swapping indices
            for j, rowNode in enumerate(tmpNodes[:-1]):
                for colNode in tmpNodes[j+1:]:
                    if colNode != mergeNode:
                        newMatrix[newNodeToIndex[rowNode], newNodeToIndex[colNode]] = distMatrix[nodeToIndex[rowNode], nodeToIndex[colNode]]
                        newMatrix[newNodeToIndex[colNode], newNodeToIndex[rowNode]] = distMatrix[nodeToIndex[colNode], nodeToIndex[rowNode]]
                    else:
                        nodeProduct = product([rowNode], digraph.adjList[colNode])
                        sumOfProductDistance = sum([distMatrix[nodeToIndex[node1], nodeToIndex[node2]]*node2.leaves for node1, node2 in nodeProduct])
                        newMatrix[newNodeToIndex[rowNode], newNodeToIndex[colNode]] = sumOfProductDistance/colNode.leaves
                        newMatrix[newNodeToIndex[colNode], newNodeToIndex[rowNode]] = sumOfProductDistance/colNode.leaves

            # Update all
            distMatrix = newMatrix
            nodeToIndex = newNodeToIndex
        
        return digraph
    
    def nj(self) -> Digraph:
        """Performs neighbor joining on distance matrix with corresponding taxa.
        """
        distMatrix = self.distMatrix
        # Initialize list of node and maps of nodes -> index
        tmpNodes = [Node(taxon) for taxon in self.taxa]
        digraph = Digraph(tmpNodes)
        nodeToIndex = {node:i for i, node in enumerate(tmpNodes)}
        
        # Find minimum
        for i in range(1, len(self.taxa)):
            # Set up neighbor joining matrix
            njMatrix = np.zeros(distMatrix.shape)
            totalDistance = distMatrix.sum(axis=1) # row sum
            for row in range(len(tmpNodes)):
                for col in range(len(tmpNodes)):
                    if row != col:
                        njMatrix[row, col] = (len(tmpNodes)-2)*distMatrix[row, col] - totalDistance[row] - totalDistance[col]

            # Find min of distMatrix and obtain their corresponding nodes and assign distance
            row, col = np.unravel_index(njMatrix.argmin(), njMatrix.shape)
            
            delta = (totalDistance[row] - totalDistance[col])/((len(tmpNodes)-2) if len(tmpNodes) > 2 else 2)
            distance = distMatrix[row, col]
            rowNode, colNode = tmpNodes[row], tmpNodes[col]
            rowNode.distance, colNode.distance = 0.5*(distance + delta), 0.5*(distance - delta)
            
            # Merge row/col Nodes into new node, remove them and append merged node
            mergeNode = digraph.join((rowNode, colNode), tmpNodes)
            print(f"Distance = {distMatrix[row, col]} +/- {delta} | Merged: {rowNode} with {colNode}") # For debugging
            
            # Reinitialize list of node and maps of nodes -> index
            newMatrix = np.zeros((len(self.taxa)-i, len(self.taxa)-i))
            newNodeToIndex = {node:i for i, node in enumerate(tmpNodes)}
            
            # Loop in upper triangle fashion, lower triangle can be filled by swapping indices
            for j, rowNode in enumerate(tmpNodes[:-1]):
                for colNode in tmpNodes[j+1:]:
                    if colNode != mergeNode:
                        newMatrix[newNodeToIndex[rowNode], newNodeToIndex[colNode]] = distMatrix[nodeToIndex[rowNode], nodeToIndex[colNode]]
                        newMatrix[newNodeToIndex[colNode], newNodeToIndex[rowNode]] = distMatrix[nodeToIndex[colNode], nodeToIndex[rowNode]]
                    else:
                        nodeProduct = product([rowNode], digraph.adjList[colNode])
                        sumOfProductDistance = sum([distMatrix[nodeToIndex[node1], nodeToIndex[node2]] for node1, node2 in nodeProduct])
                        newMatrix[newNodeToIndex[rowNode], newNodeToIndex[colNode]] = (sumOfProductDistance - distance)/2
                        newMatrix[newNodeToIndex[colNode], newNodeToIndex[rowNode]] = (sumOfProductDistance - distance)/2
            
            # Update all
            distMatrix = newMatrix
            nodeToIndex = newNodeToIndex
        
        return digraph


    