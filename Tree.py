from itertools import product

class AdjList:
    def __init__(self):
        pass
    
    def root(self) -> Node:
        pass

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
    
    def __contains__(self, other):
        pass
        
    def __hash__(self) -> int: # DO NOT HASH SOMETHING THAT IS MUTABLE
        return hash((self.taxon, self.children))
    
    def __repr__(self) -> str:
        return "-".join([leaf.taxon for leaf in self.leaves()])
    
    def setDistance(self, distance: float):
        if self.isLeaf:
            self.distance = distance
        else:
            self.distance = distance - max([node.totalDistance() for node in self.children])
        return
    
    def leaves(self) -> list:
        if self.isLeaf:
            return [self]
        leaves = []
        for node in self.children:
            leaves += node.leaves()
        return leaves
    
    def totalDistance(self) -> float:
        if self.isLeaf:
            return self.distance
        return self.distance + max([node.totalDistance() for node in self.children])
    
    def unroot(self) -> AdjList:
        pass

class Newick:
    @staticmethod
    def ToNewick(root: Node) -> str:
        """Converts rooted tree into newick string.
        """
        distance = "" if root.distance is None else f":{root.distance}"
        if root.isLeaf:
            return f"{root.taxon}{distance}"
        else:
            return f"({','.join([Newick.ToNewick(node) for node in root.children])}){distance}"
    
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
        if colonIndex == -1:
            distance = None
        else:
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
    
    @staticmethod
    def Contain(node1: Node, node2: Node) -> bool:
        """Evaluates if node2 is contained in node1, uses loose equality.
        """
        pass