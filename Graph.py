from dataclasses import dataclass

# Convert node into graph node
# Use digraph to represent a tree
# Convert all functions to use digraph

@dataclass
class Node:
    taxon: str
    distance: float = 0
    leaves: int = 1 # If node is leaf leaves = 1, else leaves = sum of children's leaves
    
    def __eq__(self, other):
        return self.taxon == other.taxon
    
    def __hash__(self):
        return hash(self.taxon)
    
    def __repr__(self):
        return f"{self.taxon}:{self.distance}"
    
class Digraph:
    def __init__(self, nodes: list):
        self.adjList = {node:[] for node in nodes}
        self._tmpNodes = nodes # private
        
    def __eq__(self, other) -> bool:
        pass
        
    def __str__(self) -> str:
        printValue = ""
        for key, item in self.adjList.items():
            nodes = ' | '.join([f"{node}" for node in item])
            printValue += f"{key}: [ {nodes} ]\n"
            
        return printValue

    def join(self, nodes: tuple, distances: tuple) -> Node:
        """Joins a tuple of nodes together
        | Distances tuple correspond to InterNode -> aNode
        | InterNode taxon are the node tuple's taxon joined with '-'
        """
        interNode = Node('-'.join([node.taxon for node in nodes]), leaves=sum([node.leaves for node in nodes]))
        self.adjList[interNode] = []
        for node, distance in zip(nodes, distances):
            node.distance = distance
            self.adjList[interNode].append(node)
            self._tmpNodes.remove(node)
        self._tmpNodes.append(interNode)
        
        return interNode
    
    def toNewick(self) -> str:
        
        pass
    
        