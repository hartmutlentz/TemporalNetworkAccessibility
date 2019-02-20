#! /usr/bin/python
#

import random
import numpy as np
try:
    import networkx as nx  # Optional, for CM model
except ImportError:
    print("Networkx-Package is not installed.\
        Configuration-Model (CM) not available.")


class TemporalEdgeList():
    """ Class for temporal edgelists as triples (u,v,d),
        Where (u,v) is an edge at time d.
    """
    def __init__(self, fname, directed, timecolumn=2):
        self.edges = np.loadtxt(fname, usecols=(0, 1, timecolumn), dtype='int')
        self.__clean_edges()
        self.is_directed = directed
        # time infos
        self.real_times = set(np.loadtxt(fname, usecols=(timecolumn,),
                                         dtype='int', unpack=True))
        self.maxtime = max(self.real_times)
        self.mintime = min(self.real_times)
        self.timespan = self.maxtime - self.mintime
        self.possible_times = range(min(self.real_times),
                                    max(self.real_times)+1)
        # init
        self.snapshots = self.__update_snapshots()
        self.static_edges = self.__get_static_edges()

    def __has_matrix_friendly_node_labels(self):
        # check if node labels are matrix friendly
        nodes1, nodes2 = zip(*self.static_edges)
        nodes = []
        nodes.extend(nodes1)
        nodes.extend(nodes2)
        nodes = set(nodes)
        nodes = list(nodes)
        nodes.sort()
        if nodes == range(len(nodes)):
            return True
        else:
            return False

    def __update_snapshots(self):
        # dict {d:[(u1,v1),(u2,v2)], ...}
        # use after self.edges have been changed/created
        t = dict([(i, []) for i in self.possible_times])

        for u, v, d in self.edges:
            t[d].append((u, v))

        return t

    def __update_edges(self):
        # reads snapshots and rewrites edgelist
        # use after self.snapshots have been changed/created
        new_edges = []
        for time in self.snapshots:
            for u, v in self.snapshots[time]:
                new_edges.append((u, v, time))
        self.edges = new_edges

    def __get_static_edges(self):
        # all edges present in the static network
        e = [(u, v) for u, v, d in self.edges]
        return list(set(e))

    def __clean_edges(self):
        # remove multiple edges
        the_edges = [(u, v, t) for (u, v, t) in self.edges]
        cut = set(the_edges)
        if len(cut) != len(the_edges):
            print('Removed multiple edges in dataset.')
        self.edges = list(cut)

    def dilute(self, p=0.5):
        """ Keep each edge in each snapshot with probability p,
            i.e. remove each such edge with probability 1-p.

        """
        assert p > 0.0, "Probability must be greater than zero."
        if p >= 1.0:
            return

        t = dict([(i, []) for i in self.possible_times])

        # get new edges
        for d in self.snapshots:
            for i, j in self.snapshots[d]:
                if random.random() < p:
                    t[d].append((i, j))

        self.snapshots = t
        self.__update_edges
        self.static_edges = self.__get_static_edges()

    def edge_occurrence_times(self):
        # dict {(u,v):[t1,t2,...],...}
        et = dict([(se, []) for se in self.static_edges])

        for u, v, d in self.edges:
            et[(u, v)].append(d)

        # sort times
        for times in et.values():
            times.sort()

        return et

    def node_occurrence_times(self):
        # dict {v:[t1, t2, ...]}
        et = dict([(se, []) for se in self.nodes()])

        for u, v, d in self.edges:
            et[u].append(d)
            et[v].append(d)

        # sort times
        for times in et.values():
            times.sort()

        return et

    def GST(self):
        # alias
        self.shuffle_snapshot_times()

    def shuffle_snapshot_times(self):
        # shuffles all snapshots
        new_keys = list((self.snapshots).keys())
        random.shuffle(new_keys)

        new_t_edges = {}
        for i in self.possible_times:
            new_t_edges[new_keys.pop()] = self.snapshots[i]

        new_edges = []
        for t in new_t_edges:
            for (u, v) in new_t_edges[t]:
                new_edges.append((u, v, t))

        self.edges = new_edges
        self.snapshots = self.__update_snapshots()

    def LST(self):
        # alias
        self.shuffle_edge_times()

    def shuffle_edge_times(self):
        """ gives every edge new occurrence times at random.
            Number of occurrences is conserved.
        """
        edge_occs = self.edge_occurrence_times()
        new_edge_occs = dict([(se, []) for se in edge_occs])

        for edge in edge_occs:
            for i in range(len(edge_occs[edge])):
                new_edge_occs[edge].append(random.randint
                                           (self.mintime, self.maxtime))

        new_edges = []
        for u, v in new_edge_occs:
            for t in new_edge_occs[(u, v)]:
                new_edges.append((u, v, t))

        self.edges = new_edges
        self.snapshots = self.__update_snapshots()

    def TR(self):
        # time reversal alias
        self.time_reversal()

    def time_reversal(self):
        # revert time stamps
        new_keys = list(range(self.mintime, self.maxtime+1))

        new_snapshots = {}
        for i in self.possible_times:
            new_snapshots[new_keys.pop()] = self.snapshots[i]

        new_edges = []
        for t in new_snapshots:
            for (u, v) in new_snapshots[t]:
                new_edges.append((u, v, t))

        self.edges = new_edges
        self.snapshots = self.__update_snapshots()

        # transpose, if network is directed
        if self.is_directed:
            for graphlet in self.snapshots:
                self.__revert_graphlet(graphlet)
            self.__update_edges()

        # update static network
        self.static_edges = self.__get_static_edges()

    def __revert_graphlet(self, time):
        # reverts/transposes a single snapshot
        edges = self.snapshots[time][:]
        new_edges = [(v, u) for (u, v) in edges]
        self.snapshots[time] = new_edges

    def __graphlet_configuration_model(self, G_in):
        """ Converts network input network into graph sequence and
            generates new configuration graph.
            Number of edges is not conserved!
        """
        if G_in.is_directed():
            # inseq = list(G_in.in_degree().values())
            inseq = [k for k, v in G_in.in_degree]
            # outseq = list(G_in.out_degree().values())
            outseq = [k for k, v in G_in.out_degree]

            H = nx.directed_configuration_model(inseq, outseq)
            H = nx.DiGraph(H)
            H.remove_edges_from(H.selfloop_edges())

        else:
            seq = G_in.degree().values()

            H = nx.configuration_model(seq)
            H = nx.Graph(H)
            H.remove_edges_from(H.selfloop_edges())

        return H.edges()

    def __graphlet_to_nx_graph(self, time):
        """ Converts a snapshot to nx.(Di)Graph """
        if self.is_directed:
            G = nx.DiGraph()
        else:
            G = nx.Graph()
        G.add_nodes_from(self.nodes())
        G.add_edges_from(self.snapshots[time])
        return G

    def __randomize_graphlet(self, time, maxiterations=100):
        """
            Returns a randomized version of a graph or digraph.
            The degree sequence is conserved.
        """
        iterations = len(self.snapshots[time]) * maxiterations

        def get_legal_edgepair(ed):
            # returns a disjoint pair of edges
            def are_disjoint(fi, se):
                # condition for disjoint edges
                if fi[0] == fi[1] or se[0] == se[1]\
                    or fi[0] == se[0] or fi[0] == se[1]\
                        or fi[1] == se[0] or fi[1] == se[1]:
                            return False
                else:
                    return True

            def pair_not_in_G(x, y, eds):
                # True, if exchanged edges are not already in G.
                if (x[0], y[1]) not in eds and (y[0], x[1]) not in eds:
                    return True
                else:
                    return False

            for i in range(100*len(ed)):
                first = random.sample(ed, 1)[0]  # random.choice(ed)
                second = random.sample(ed, 1)[0]  # random.choice(ed)
                if are_disjoint(first, second) and \
                        pair_not_in_G(first, second, ed):
                        return (first, second)
            return False

        def legal_graph_condition(e):
            # conditions for useful edgelists
            if len(e) < 2:
                return False

            nodeset = set()
            for u, v in e:
                nodeset.add(u)
                nodeset.add(v)
                if len(nodeset) > 3:
                    return True

        # switch edges
        edges = set(self.snapshots[time][:])

        if legal_graph_condition(edges):
            for i in range(iterations):
                erfolg = get_legal_edgepair(edges)
                if erfolg:
                    x, y = erfolg[0], erfolg[1]

                    edges.remove((x[0], x[1]))
                    edges.remove((y[0], y[1]))

                    edges.add((x[0], y[1]))
                    edges.add((y[0], x[1]))
                # print('remaining: ', iterations-i)

        self.snapshots[time] = list(edges)

    def configuration_model(self):
        # alias
        self.CM()

    def CM(self):
        """ Configuration model.
            Faster than RE, but number of edges is not conserved.
            Use RE for smaller networks only.
        """
        for i, t in enumerate(self.snapshots):
            print("Configuration model for t= ", i, " of ", self.timespan)
            new_edges = self.\
                __graphlet_configuration_model(self.__graphlet_to_nx_graph(t))
            self.snapshots[t] = new_edges

        self.__update_edges()
        self.static_edges = self.__get_static_edges()

    def RE(self, maxiterations=100):
        # alias
        self.randomize_edges(maxiterations)

    def randomize_edges(self, maxiterations=100):
        """ Edge randomization for each graphlet
        """
        for i, j in enumerate(self.snapshots):  # range(self.maxtime):
            print("Randomizing ", i, " of ", self.timespan)
            self.__randomize_graphlet(j, maxiterations)
        self.__update_edges()
        self.static_edges = self.__get_static_edges()

    def number_of_nodes(self):
        return len(self.nodes())

    def nodes(self):
        # the nodes of the network as set.
        nodes = []
        for (u, v) in self.static_edges:
            nodes.append(u)
            nodes.append(v)
        nodes = set(nodes)
        return nodes

    def write(self, fname):
        """ writes self to txtfile.
        """
        arr = self.edges
        g = open(fname, 'w+')
        for i in range(len(arr)):
            wstring = ''
            for j in range(1, len(arr[i])):
                wstring += '\t' + str(arr[i][j])
            g.writelines((str(arr[i][0]) + wstring + '\n'))
        g.close()

        return

    def average_size(self):
        """ average edge density """
        all_edges = 0
        for time in self.snapshots:
            all_edges += len(self.snapshots[time])

        return float(all_edges) / (self.maxtime - self.mintime)

    def random_times_uniform(self):
        """ Times at random from uniform distribution
            SLOW!
        """
        prob = self.average_size() / len(self.static_edges)
        # print(prob)
        for i in range(self.mintime, self.maxtime):
            print("Random times uniform. Step ", i, " of ", self.maxtime)
            edges = []
            for e in self.static_edges:
                if random.random() < prob:
                    edges.append(e)
            self.snapshots[i] = edges
        self.__update_edges()

    def RT(self):
        # alias
        self.random_times()

    def random_times(self):
        """ Keeps the distribution of graph sizes.
            Example:
            before: time_1: 2 edges, time_2: 5 edges, time_3: 4 edges
            edges are elements of static edges (fixed):
            after: time_1: 5 edges, time_2: 4 edges, time_3: 2 edges
            edges are chosen randomly from static graph
        """
        sizes = dict([(i, len(self.snapshots[i])) for i in self.snapshots])
        # new permutation of edge densities
        new_keys = list(sizes.keys())
        random.shuffle(new_keys)

        new_sizes = {}
        for i in sizes:
            new_sizes[new_keys.pop()] = sizes[i]

        timespan = self.maxtime - self.mintime
        for i, j in enumerate(sizes):
            print("Random times. Step ", i, " of ", timespan)
            edges = set()
            # while len(edges)<new_sizes[j]:
            for _ in range(len(self.static_edges)):
                if len(edges) >= new_sizes[j]:
                    break
                edges.add(random.choice(self.static_edges))

            self.snapshots[j] = list(edges)
        self.__update_edges()


if __name__ == "__main__":
    from pprint import pprint
    the_file = '../edgelists/Test.dat'
    E = TemporalEdgeList(the_file, True, timecolumn=2)

    print(E.snapshots[0])
    E.dilute(0.09)
    print(E.snapshots[0])
    #print(len(E.edges))
    #E=TemporalEdgeList("sociopatterns_113.dat",False)
    #pprint(E.edges)
    #E.randomize_edges()
    #E.random_times()
    #print(E.average_size(),len(E.snapshots))
    #print(len(E.edges))

    #E.write("out1_RE.txt")

    #print(E.edge_occurrence_times())
    #print(E.shuffle_edge_times(E.edge_occurrence_times()))
