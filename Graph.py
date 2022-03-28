from dataclasses import dataclass

# Convert node into graph node
# Use digraph to represent a tree
# Convert all functions to use digraph

@dataclass
class Node:
    taxon: str
    distance: float = 0.0
    total: float = 0.0
    leaves: int = 1 # If node is leaf leaves = 1, else leaves = sum of children's leaves
    
    def __eq__(self, other):
        return self.taxon == other.taxon
    
    def __len__(self):
        return len(self.taxon)
    
    def __hash__(self):
        return hash(self.taxon)
    
    def __repr__(self):
        return f"({self.taxon}:{self.distance})"
    
class Digraph:
    def __init__(self, nodes: list, edges: list=None):
        self.adjList = {node:[] for node in nodes}
        
    def __eq__(self, other) -> bool:
        
        pass
        
    def __str__(self) -> str:
        printValue = ""
        for key, item in self.adjList.items():
            nodes = ' | '.join([f"{node}" for node in item])
            printValue += f"{key}: [ {nodes} ]\n"
            
        return printValue

    def join(self, nodes: tuple, tmpNodes: list) -> Node:
        """Joins a tuple of nodes together
        """
        interName = ','.join([f"{node.taxon}:{node.distance}" for node in nodes])
        interNode = Node(f"({interName})", leaves=sum([node.leaves for node in nodes]))
        self.adjList[interNode] = []
        for node in nodes:
            self.adjList[interNode].append(node)
            tmpNodes.remove(node)
        tmpNodes.append(interNode)
        
        return interNode
    
    def toNewick(self) -> str:
        nodes = list(self.adjList.keys())
        newickNode = max(nodes, key=len)
        print(f"{newickNode.taxon}:{newickNode.distance}")
        return f"{newickNode.taxon}:{newickNode.distance}"
    
        