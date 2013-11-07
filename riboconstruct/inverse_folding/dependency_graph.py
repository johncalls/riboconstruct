import collections
import random


class DependencyGraph(object):
    """ DependencyGraph is used to generate dependency subgraphs. """

    def __init__(self):
        self.connected_nodes = collections.defaultdict(set)

    def add_edge(self, edge):
        node_a, node_b = edge
        # check the number of known outgoing edges for each new node
        if len(self.connected_nodes[node_a]) == 2:
            raise ValueError("Pos %i has already two outgoing edges." % node_a)
        if len(self.connected_nodes[node_b]) == 2:
            raise ValueError("Pos %i has already two outgoing edges." % node_b)
        # check if the same edge has been added previously
        if node_a in self.connected_nodes[node_b]:
            return
        if node_b in self.connected_nodes[node_a]:
            return
        self.connected_nodes[node_a].add(node_b)
        self.connected_nodes[node_b].add(node_a)

    def iter_subgraphs(self):
        def circular_check(node, node_start, subgraph):
            """
            Checks whether a given subgraph is circular by traversing all edges
            defining the subgraph. If at some point it ends up at the start node
            the subgraph is circular otherwise it is linear. All traversed nodes
            in the subgraph are added to 'subgraph'.
            """

            subgraph.add_edge((node, node_start))
            node_prev = node_start
            while node != node_start:
                try:
                    node_a, node_b = self.connected_nodes[node]
                except ValueError:
                    # only one outgoing edge, i.e. the subgraph is acyclic
                    return False
                finally:
                    del self.connected_nodes[node]
                if node_a != node_prev:
                    node_prev = node
                    node = node_a
                else:
                    node_prev = node
                    node = node_b
                subgraph.add_edge((node_prev, node))
            return True

        while len(self.connected_nodes):
            node, connected_nodes = self.connected_nodes.popitem()

            subgraph = DependencySubgraph()
            # if there is only one outgoing edge the subgraph is acyclic
            if len(connected_nodes) == 1:
                subgraph.is_cyclic = False
                # use circular check to traverse the subgraph and add the
                # nodeitions to 'subgraph'
                circular_check(connected_nodes.pop(), node, subgraph)
                yield subgraph
            else:
                node_a, node_b = connected_nodes
                # check if the graph is circular: by doing this all traversed
                # nodeitions are added to 'subgraph'
                if circular_check(node_a, node, subgraph):
                    subgraph.is_cyclic = True
                    yield subgraph
                else:
                    subgraph.is_cyclic = False
                    # if the subgraph is acyclic the nodeitions of the other
                    # 'leg' of the subgraph have to be added as well
                    circular_check(node_b, node, subgraph)
                    yield subgraph


class DependencySubgraph(object):
    """
    A basepair (i,j) of a RNA structures is dependent on the basepair (k,l) of
    another structure iff either i==l or j==k. A dependency subgraph is a
    concatination of dependent basepair nodes. Subgraphs can either be
    circular (e.g. (i-j-k-l-i)) or linear (e.g. (i-j-k-l)).

    DependencySubgraph provides diverse methods to traverse such subgraphs.
    """

    def __init__(self):
        self.connected_nodes = collections.defaultdict(list)

    def add_edge(self, edge):
        node_a, node_b = edge
        # check the number of known outgoing edges for each new node
        if len(self.connected_nodes[node_a]) == 2:
            raise ValueError("Pos %i has already two outgoing edges." % node_a)
        if len(self.connected_nodes[node_b]) == 2:
            raise ValueError("Pos %i has already two outgoing edges." % node_b)
        self.connected_nodes[node_a].append(node_b)
        self.connected_nodes[node_b].append(node_a)

    def iter_nodes(self):
        return self.connected_nodes.iterkeys()

    def iter_nodes_sorted(self):
        if self.is_cyclic:
            for node in self.iter_nodes():
                yield node
        else:
            nodes_iter = self.iter_nodes()
            node_start = nodes_iter.next()
            while len(self.connected_nodes[node_start]) == 2:
                node_start = nodes_iter.next()
            yield node_start
            node_prev = node_start
            node = self.connected_nodes[node_start][0]
            while len(self.connected_nodes[node]) != 1:
                yield node
                node_a, node_b = self.connected_nodes[node]
                if node_a != node_prev:
                    node_prev = node
                    node = node_a
                else:
                    node_prev = node
                    node = node_b
            yield node

    def iter_edges(self):
        nodes_iter = self.iter_nodes()
        node_start = nodes_iter.next()
        if not self.is_cyclic:
            while len(self.connected_nodes[node_start]) == 2:
                node_start = nodes_iter.next()
        node_prev = node_start
        node = self.connected_nodes[node_start][0]
        while node != node_start and len(self.connected_nodes[node]) != 1:
            yield (node_prev, node)
            node_a, node_b = self.connected_nodes[node]
            if node_a != node_prev:
                node_prev = node
                node = node_a
            else:
                node_prev = node
                node = node_b
        yield (node_prev, node)

    def get_nodes_iter(self, node_start, rand=False):
        def iter_circular_subgraph(node_start):
            yield node_start
            node_prev = node_start
            # randomly choose in which direction to traverse the subgraph
            if rand:
                random.shuffle(self.connected_nodes[node_start])
            node_a, node_b = self.connected_nodes[node_start]
            if node_a != node_prev:
                node = node_a
            else:
                node = node_b
            yield node
            # traverse the subgraph
            while node != node_start:
                node_a, node_b = self.connected_nodes[node]
                if node_a != node_prev:
                    node_prev = node
                    node = node_a
                else:
                    node_prev = node
                    node = node_b
                yield node

        def iter_linear_subgraph(node_start, node_prev):
            yield node_prev
            yield node_start
            node = node_start
            # traverse the 'leg'
            while True:
                try:
                    node_a, node_b = self.connected_nodes[node]
                except ValueError:
                    return
                if node_a != node_prev:
                    node_prev = node
                    node = node_a
                else:
                    node_prev = node
                    node = node_b
                yield node

        if self.is_cyclic:
            return iter_circular_subgraph(node_start)
        # if the subgraph is acyclic then return either the full 'leg' (if the
        # start node is at either end of the subgraph) or both 'legs' (with the
        # the start node as root)
        else:
            if len(self.connected_nodes[node_start]) == 2:
                # two outgoing edges available; randomly select which one to
                # traverse first
                if rand:
                    random.shuffle(self.connected_nodes[node_start])
                node_a, node_b = self.connected_nodes[node_start]
                return (iter_linear_subgraph(node_a, node_start),
                        iter_linear_subgraph(node_b, node_start))
            else:
                # if only one outgoing edge available
                node_a = self.connected_nodes[node_start][0]
                return (iter_linear_subgraph(node_a, node_start),)

    def get_rand_nodes_iter(self):
        i = random.randint(0, len(self.connected_nodes) - 1)
        node_start = list(self.iter_nodes())[i]
        return self.get_nodes_iter(node_start, rand=True)
