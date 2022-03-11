import numpy as np
from itertools import product

class Node:
    def __init__(self, taxon: str=None, children: tuple=(), distance: float=0.0):
        self.taxon = taxon
        self.children = children
        self.distance = distance
        if self.taxon is not None:
            self.isLeaf = True
        else:
            self.isLeaf = False
            
    def __len__(self) -> int:
        if self.isLeaf:
            return 1
        length = 0
        for node in self.children:
            length += len(node)
        return length
    
    def __eq__(self, other, strict=True):
        """Equality of two nodes |
        If strict, check equal distance |
        else loose, distance is not checked.
        """
        equalDistance = (self.distance == other.distance) if strict else True
        
        if self.isLeaf ^ other.isLeaf:
            return False
        if self.isLeaf and other.isLeaf:
            return self.taxon == other.taxon and equalDistance
        
        flag = {node:False for node in self.children}
        for (selfChild, otherChild) in product(self.children, other.children):
            flag[selfChild] = flag[selfChild] or selfChild.__eq__(otherChild, strict)
                
        return all(flag.values()) and equalDistance
        
    def __hash__(self) -> int: # DO NOT HASH SOMETHING THAT IS MUTABLE
        return hash((self.taxon, self.children))
    
    def __repr__(self) -> str:
        return "-".join(self.leaves())
    
    def setDistance(self, distance: float):
        if self.isLeaf:
            self.distance = distance
        else:
            self.distance = distance - max([node.totalDistance() for node in self.children])
        return
    
    def leaves(self) -> list:
        if self.isLeaf:
            return [self.taxon]
        leaves = []
        for node in self.children:
            leaves += node.leaves()
        return leaves
    
    def totalDistance(self) -> float:
        if self.isLeaf:
            return self.distance
        return self.distance + max([node.totalDistance() for node in self.children])

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

class Newick:
    @staticmethod
    def ToNewick(root: Node) -> str:
        """Converts rooted tree into newick string.
        """
        if root.isLeaf:
            return f"{root.taxon}:{root.distance}"
        else:
            return f"({','.join([Newick.ToNewick(node) for node in root.children])}):{root.distance}"
    
    @staticmethod
    def ToTree(newick: str) -> Node:
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

        colonIndex = newick.rfind(':')
        distance = float(newick[colonIndex+1:])
        newick = newick[:colonIndex]
        if newick[0] != '(':
            return Node(taxon=newick, distance=distance)
        
        commaList = findComma(newick[:colonIndex])
        newickList = []
        initial = 1
        for commaIndex in commaList:
            newickList.append(newick[initial:commaIndex])
            initial = commaIndex+1
        
        nodeList = []
        for newick in newickList:
            nodeList.append(Newick.ToTree(newick))

        return Node(children=tuple(nodeList), distance=distance)
    
    @staticmethod
    def Equal(node1: Node, node2: Node, strict=True) -> bool:
        """Equality of two nodes |
        If strict, check equal distance |
        else loose, distance is not checked.
        """
        return node1.__eq__(node2, strict)
    
    @staticmethod
    def Clade(node1: Node, node2: Node) -> bool:
        """Evaluates if node2 is a clade of node1, uses loose equality.
        """
        if Newick.Equal(node1, node2, strict=False):
            return True
        else:
            flag = False
            for child in node1.children:
                flag = flag or Newick.Clade(child, node2)
        return flag
