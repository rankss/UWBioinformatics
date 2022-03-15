import numpy as np
from itertools import product
from Tree import Node

class Cluster:
    def __init__(self, distMatrix: np.ndarray, taxa: list):
        self.distMatrix = np.array(distMatrix, dtype=float)
        np.fill_diagonal(self.distMatrix, np.inf)
        self.taxa = taxa
        
    def upgma(self) -> Node:
        """Performs UPGMA on a distance matrix with corresponding taxa.
        """
        # https://codereview.stackexchange.com/questions/263416/upgma-tree-building-in-python with a modifications
        distMatrix = self.distMatrix
        # Initialize list of node and maps of nodes -> index
        nodes = [Node(taxon=taxon) for taxon in self.taxa]
        nodeToIndex = {node:i for i, node in enumerate(nodes)}
        
        for i in range(1, len(self.taxa)):
            # Find min of distMatrix and obtain their corresponding nodes and assign distance
            index = np.unravel_index(distMatrix.argmin(), distMatrix.shape)
            distance = distMatrix[index[0], index[1]]
            rowNode, colNode = nodes[index[0]], nodes[index[1]]
            rowNode.setDistance(distance/2), colNode.setDistance(distance/2)
            # Merge row/col Nodes into new node, remove them and append merged node
            mergeNode = Node(children=(rowNode, colNode))
            nodes.remove(rowNode), nodes.remove(colNode)
            nodes.append(mergeNode)
            
            # Reinitialize list of node and maps of nodes -> index
            newMatrix = np.full((len(self.taxa)-i, len(self.taxa)-i), np.inf)
            newNodeToIndex = {node:i for i, node in enumerate(nodes)}
                
            # Loop in upper triangle fashion, lower triangle can be filled by swapping indices
            for j, rowNode in enumerate(nodes[:-1]):
                for colNode in nodes[j+1:]:
                    if colNode != mergeNode:
                        newMatrix[newNodeToIndex[rowNode], newNodeToIndex[colNode]] = distMatrix[nodeToIndex[rowNode], nodeToIndex[colNode]]
                        newMatrix[newNodeToIndex[colNode], newNodeToIndex[rowNode]] = distMatrix[nodeToIndex[colNode], nodeToIndex[rowNode]]
                    else:
                        nodeProduct = product([rowNode], list(colNode.children))
                        result = np.array([[distMatrix[nodeToIndex[rowNode], nodeToIndex[colNode]]*len(colNode), len(colNode)] for rowNode, colNode in nodeProduct])
                        newMatrix[newNodeToIndex[rowNode], newNodeToIndex[colNode]] = result[:,0].sum()/result[:,1].sum()
                        newMatrix[newNodeToIndex[colNode], newNodeToIndex[rowNode]] = result[:,0].sum()/result[:,1].sum()

            # Update all
            distMatrix = newMatrix
            nodeToIndex = newNodeToIndex
        
        return nodes[0]
    
    def neighborJoining(self) -> Node:
        # Charlotte TODO
        
        pass


    