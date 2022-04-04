from dataclasses import dataclass

@dataclass
class Node:
    taxon: str
    distance: float = 0.0
    total: float = 0.0
    leaves: int = 1 # If node is leaf leaves = 1, else leaves = sum of children's leaves
    
    def __eq__(self, other) -> bool:
        return self.taxon == other.taxon and self.distance == other.distance
    
    def __len__(self) -> int:
        return len(self.taxon)
    
    def __hash__(self) -> int:
        return hash(self.taxon)
    
    def __repr__(self) -> str:
        return f"({self.taxon}:{self.distance})"

@dataclass
class Edge:
    node1: Node
    node2: Node
    weight: float
    
class Digraph:
    def __init__(self, nodes: list[Node], edges: list=None):
        self.adjList = {node:[] for node in nodes}
        if edges:
            for edge in edges:
                edge.node2.distance = edge.weight
                self.adjList[edge.node1].append()
        
    def __eq__(self, other) -> bool:
        for key in self.adjList:
            if not (key in other.adjList and set(self.adjList[key]) == set(other.adjList[key])):
                return False
            
        return True
        
    def __str__(self) -> str:
        printValue = ""
        for key, item in self.adjList.items():
            nodes = ' | '.join([f"{node}" for node in item])
            printValue += f"{key}: [ {nodes} ]\n"
            
        return printValue

    def join(self, nodes: tuple[Node], tmpNodes: list[Node]) -> Node:
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
        return f"{newickNode.taxon}:{newickNode.distance}"
    