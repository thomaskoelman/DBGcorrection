

class Graph:
    def __init__(self):
        self.nodes = list()
        self.edges = list()

    def has_node(self, node):
        return node in self.nodes

    def add_node(self, node):
        if not self.has_node(node):
            self.nodes.append(node)
        return self

    def has_edge(self, tuple):
        return tuple in self.edges

    def add_edge(self, tuple):
        if self.has_node(tuple[0]) and self.has_node(tuple[1]) and not self.has_edge(tuple):
            self.edges.append(tuple)
        return self

    def delete_node(self, node):
        if node in self.nodes:
            self.nodes.remove(node)
            edges_out = list(filter(lambda edge: (edge[0] == node), self.edges))
            edges_in = list(filter(lambda edge: edge[1] == node, self.edges))
            for edge in edges_out:
                if edge in self.edges:
                    self.edges.remove(edge)
            for edge in edges_in:
                if edge in self.edges:
                    self.edges.remove(edge)

    def delete_edge(self, tuple):
        if tuple in self.edges:
            self.edges.remove(tuple)
        return self

    def count_edges(self, node, incoming):
        i = 1 if incoming else 0
        edges = list(filter(lambda edge: edge[i] == node, self.edges))
        return len(edges)

    def get_neighbours(self, node, incoming=False):
        i = 1 if incoming else 0
        neighbours = list()
        edges = list(filter(lambda edge: edge[i] == node, self.edges))
        for edge in edges:
            neighbours.append(edge[1-i])
        return neighbours

    def get_simple_paths(self):
        def get_rest(node):
            is_branching = self.count_edges(node, incoming=False) > 1
            is_before_merge = bool([n for n in self.get_neighbours(node, incoming=False) if self.count_edges(n, incoming=True) > 1])
            is_ending = self.count_edges(node, incoming=False) == 0
            if is_branching or is_ending or is_before_merge:
                return list()
            else:
                neighbours = self.get_neighbours(node)
                next = neighbours[0]
                return list([next]) + get_rest(next)

        paths = list()
        c = 0
        l = len(self.nodes)
        for node in self.nodes:
            is_merging = self.count_edges(node, incoming=True) > 1
            is_starting = self.count_edges(node, incoming=True) == 0
            is_after_branch = bool([n for n in self.get_neighbours(node, incoming=True) if self.count_edges(n, incoming=False) > 1])
            one_edge_out = self.count_edges(node, incoming=False) == 1
            if (is_starting and one_edge_out) or (is_merging and one_edge_out) or (is_after_branch and one_edge_out):
                #what if there is a perfect cycle?
                path = list([node]) + get_rest(node)
                if len(path) > 1:
                    paths.append(path)
        return paths


# g = Graph()
# g.add_node('ABCD')
# print(g.has_node('ABCD'))
# g.add_node('DCBA')
# print(g.has_node('DCBA'))
# g.add_edge(('ABCD','DCBA'))
# print(g.has_edge(('ABCD','DCBA')))
# print(g.has_edge(('ABCD','CDBA')))
# g.add_node('E')
# g.add_edge(('ABCD','E'))
# print(g.count_edges('ABCD', incoming=False))
# print(g.has_edge(('ABCD','E')))
# print(g.delete_edge(('ABCD','E')))
# print(g.has_edge(('ABCD','E')))
# g.add_node('F')
# g.add_edge(('E','F'))
# print(g.get_simple_paths())
# g.add_edge(('ABCD','E'))
# print(g.get_simple_paths())
# g.add_node('H')
# g.add_edge(('F','H'))
# print(g.get_simple_paths())
# g.delete_node('DCBA')
# print(g.get_simple_paths())
