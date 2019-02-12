#! /usr/bin/python

from scipy.sparse import csr_matrix, lil_matrix
import scipy.sparse as sp
import scipy.stats
import numpy as np
from numpy import loadtxt, zeros
import scipy
from collections import defaultdict
from scipy.io import mmread
import random
import copy
import itertools


class AdjMatrixSequence(list):
    """
    list of sparse matrices.

    class inherits from list.
    The constructor expects filename of u,v,w,d-edgelist.

    It constructs a list, where each entry is
    a sparse matrix (default: csr).
    Matrices in the list are ordered by time.

    The indices of the matrix represent the nodes. Nodes are reindexed to
    [0..number_of_nodes].

    """
    def __init__(self, edgelist_fname, directed, write_label_file=False,
                 columns=(0, 1, 2), firsttime=None, lasttime=None):
        list.__init__(self)
        self.first_day = firsttime
        self.last_day = lasttime
        self.fname = edgelist_fname
        self.cols = columns
        self.label_file = write_label_file
        self.is_directed = directed

        self.matricesCreation()
        if not self.is_directed:
            self.as_undirected()
        self.number_of_nodes = scipy.shape(self[0])[0]

    def copy(self):
        """ alias """
        return copy.copy(self)

    def deepcopy(self):
        """ alias """
        return copy.deepcopy(self)

    def __scipy_version_for_large_matrices(self, ref="0.14.0"):
        """
        Checks if the installed version of Scipy is at least ref.
        This is necessary for handling very large sparse matrices.
        For scipy versions < 0.14.0. the variables indptr and indices are
        stored as int32. Thus, the number of non-zero entries is restricted to
        2^31 \sim 10^9 elements.
        Returns True, if installed Version > ref. False otherwise.

        """
        # split versions into numbers
        x = scipy.__version__.split(".")
        y = ref.split(".")

        # compare adjusted lengths
        if len(x) == len(y):
            return x >= y
        elif len(x) < len(y):
            return x >= y[:len(x)]
        else:
            return x[:len(y)] >= y

    def info_scipy_version(self):
        """ Print information about scipy version and maximum Matrix size. """
        if self.__scipy_version_for_large_matrices():
            print("Scipy version can handle matrices with more than " +
                   "2^31 nonzero elements.")
        else:
            print("Number of nonzero elements is restricted to 2^31.")

    def groupByTime(self, edges):
        """ returns list of tupels: [(d,[(u,v),...]),...]. """
        dct = defaultdict(list)
        for u, v, d in edges:
            dct[d].append((u, v))
        dct_s = dict.fromkeys(range(0, self.last_day-self.first_day), [])
        for d in dct:
            dct_s[d-self.first_day] = dct[d]
        return dct_s.items()

    def reindex(self, edges):
        """ for any index in edges returns dict with new_indices
            [0...number_of_unique_indices].
        """
        us, vs, ds = zip(*edges)
        nodes = set(us) | set(vs)
        old_to_new = {}
        for i, n in enumerate(nodes):
            old_to_new[n] = i
        return old_to_new

    def all_time_windows(self, max_window=None):
        """ summation over all time windows """
        p_corr = {}
        if not max_window:
            mw = len(self)

        for twindow in range(mw):
            print(twindow)
            p_corr[twindow] = self.single_time_window(twindow)
        return p_corr

    def single_time_window(self, windowsize):
        """ summation over a time window """
        prop_corr = []

        for i in range(len(self) - windowsize):
            prop_corr.append(self.two_link_density(i, i+windowsize))

        return (scipy.stats.scoreatpercentile(prop_corr, 25),
                scipy.stats.scoreatpercentile(prop_corr, 50),
                scipy.stats.scoreatpercentile(prop_corr, 75))

    def two_link_density(self, index1, index2, norm=True):
        """ the link density for 2 step paths """
        C = self[index1] * self[index2]

        nlinks = C.sum()
        n = self.number_of_nodes

        if norm:
            return float(nlinks) / (float((n - 2) * (n**2 - n)) +
                                    float(n**2 - n))
        else:
            return float(nlinks) / float((n**2 - n))

    def deep_product(self, twindow=1, start=0):
        """ Product A_1*A_7*A_14... """
        C = self[start].copy()
        links = {}
        links[start] = (C * C).sum()

        for i in range(start+twindow, len(self)-twindow, twindow):
            C = C * self[i]
            links[i] = C.nnz()
            if C.nnz == 0:
                break

        return links

    def __matrix_mean_degree(self, A):
        """ the mean degree of A. """
        return float(A.nnz) / (self.number_of_nodes - 1)

    def __matrix_LCC_size(self, A):
        """ The size of the L(S)CC of a matrix.
            Note that the matrices here are interpreted as (bi)directed nets.
        """
        n, l = sp.csgraph.connected_components(A, connection='strong')
        lcc_size = np.max(np.bincount(l))
        if lcc_size == 1:
            return 0.0
        else:
            return lcc_size / float(self.number_of_nodes)

    def average_path_length(self, M, diameter):
        """ returns average pathlength of a snapshot, where all path
            lengths > 1 are considered.
        """
        x = [(M**i).nnz for i in range(2, diameter+1)]

        return np.mean(x)

    def long_paths_per_snapshot(self, max_path_length):
        """ Computes the number of paths longer than 1 in each snapshot.
            This issue has been discussed in
                Grindrod et al.
                *Communicability across evolving networks*
                Phys. Rev. E, 2011.

            Matrices in At should be corrected, if many entries are > 1.
            max_path_length should be the diameter of the aggregated network.
        """
        d = {}
        for i in range(len(self)):
            print(i)
            p = self.average_path_length(self[i], max_path_length)
            # if p > 1:
            d[i] = p

        return d

    def long_path_correction(self, diameter):
        """ replace A_i by \sum _{i=1} ^diameter A_i. This takes into account
            paths of length > 1 in each snapshot.
            See
                Grindrod et al.
                *Communicability across evolving networks*
                Phys. Rev. E, 2011.
        """
        for i in range(len(self)):
            print("Correcting snapshots for long paths. Step ", i)
            M = self[i].copy()

            for j in range(2, diameter):
                M = M + self[i]**j

            self[i] = M

        print("---> paths up to length ", diameter, \
            " are now considered in snapshots.")

        return

    def LCCs(self):
        """ returns information about the size of the LCC and the mean degree
            of the snapshots.
        """
        crit = {}
        for i, M in enumerate(self):
            print("LCC ", i)
            crit[i] = (self.__matrix_mean_degree(M), self.__matrix_LCC_size(M))

        return crit

#    def path_density_of_A(self, A):
#        """ The path density of an Adjacency Matrix A """
#        paths = 0
#        n = scipy.shape(A)[0]
#
#        for i in range(n):
#            out_size =\
#                len(sp.csgraph.depth_first_order(A,
#                    i, return_predecessors=False))-1
#            paths += out_size
#
#        return float(paths) / n**2
    def LSCC_nodes(self, A):
        """ Return nodes of the largest stringly connected component of a
            network given by a sparse matrix A.

            Parameters
            ----------
            A - scipy.sparse matrix
                Adjacency matrix of the network. Assumed to be directed.

            Returns
            -------
            list
                Node labels of the nodes in the LSCC.
        """
        n, labs = sp.csgraph.connected_components(A, connection="strong")

        nodes = set()
        for i in range(len(labs)):
            if labs[i] == 0:
                nodes.add(i)

        return nodes

    def path_density_of_A(self, A, normalize=False):
        """ Return the path density of a network given by adjacency matrix A.
            Selfloops are explicitely included.
        """
        n = scipy.shape(A)[0]
        all_nodes = set(range(n))

        lscc_nodes = self.LSCC_nodes(A)
        rest_nodes = all_nodes - lscc_nodes

        paths = 0

        # LSCC nodes
        lscc_sample_node = next(iter(lscc_nodes))
        lscc_range = len(sp.csgraph.depth_first_order(A,
                         lscc_sample_node, return_predecessors=False))
        for i in lscc_nodes:
            paths += lscc_range

        # Remaining nodes
        for i in rest_nodes:
            rng = len(sp.csgraph.depth_first_order(A,
                      i, return_predecessors=False))
            paths += rng

        if normalize:
            return float(paths) / float(n**2)
        else:
            return paths

    def static_path_density(self, normalize=False):
        """
        Path density of the aggregated network. This algorithm first
        computes the range of every node, i.e. the size of its
        out-component. The path density is then given by the sum over all
        node ranges.

        The method makes use of the giant component structure of the
        network to save cpomutation time.

        An explanation is given in:
        Lentz, H. H. K. et al.
        *Disease Spread through Animal Movements: A Static and Temporal
        Network Analysis of Pig Trade in Germany.*
        PLOS ONE 11, e0155196-32 (2016).

        Parameters
        ----------
        normalize - Boolean
            If True, the path density is given by the number of all paths
            normalized by N**2, where N is the number of nodes.

        Returns
        -------
        int or float.
            The number of paths in the static network. If normalized, the
            path density.

        Example
        -------
        >>> # At = AdjMatrixSequence(<...>)
        >>> # This method actually works in an At-Object only.
        >>> A = sp.lil_matrix((5,5))
        >>> A[0,1] = 1
        >>> A[4,0] = 1
        >>> A[1,2] = 1
        >>> A[2,3] = 1
        >>> A[3,1] = 1
        >>> print static_path_density(A, normalized=False)
            18

        Note that the network is considered directed. There are 18 paths,
        i.e. 13 paths between nodes and 5 paths (selfloops) are given by
        definition.
        """
        return self.path_density_of_A(self.cumulated(), normalize)

    def step_by_step_static_path_density(self, ende=None):
        """ Returns list. [index=Aggregation depth: static path density]

        """
        pd = []
        for i, Cn in enumerate(self.step_by_step_aggregation(ende)):
            print('Static path density. Step ', i)
            pd.append(self.path_density_of_A(Cn))

        return pd

    def step_by_step_aggregation(self, ende=None):
        """ Returns matrix list of all aggregated networks,
            i.e. [A1, A1+A2, A1+A2+A3, ...]

        """
        C = csr_matrix((self.number_of_nodes, self.number_of_nodes),
                       dtype=np.int32)
        li = []

        if ende:
            e = ende
        else:
            e = len(self)

        for i in range(e):
            C = C+self[i]
            li.append(C)

        return li

    def cumulated(self, weighted=False, start=0, ende=None):
        """ Returns Cumulated Graph as Matrix """
        C = csr_matrix((self.number_of_nodes, self.number_of_nodes),
                       dtype=np.int32)

        if ende:
            e = ende
        else:
            e = len(self)

        for i in range(start, e):
            C = C + self[i]

        if weighted:
            return C
        else:
            self.bool_int_matrix(C)
            return C

    def coarse_grain(self, aggregate, return_copy=False):
        """ coarse grain the list, i.e. partial aggregation of the network.
            aggregate - gives the number of matrices to be summed up.
            In numbers to get the idea:

                coarse_grain([1,2,3,4,5,6], 2) = [1+2, 3+4, 5+6].
             or coarse_grain([1,2,3,4,5,6], 3) = [1+2+3, 4+5+6].

            Finite Size effects are ignored! Thus
            coarse_grain([1,2,3,4,5,6], 2) gives the same result as
            coarse_grain([1,2,3,4,5,6,7], 2), since the last element cannot
            be aggregated into a 2-aggregate.
        """
        def main_loop(x):
            new_list = []
            while x:
                partial = []
                for i in range(aggregate):
                    try:
                        partial.append(x.pop(0))
                    except IndexError:
                        return new_list
                new_list.append(sum(partial))
            return new_list

        assert aggregate <= len(self),\
            'Aggregate must contain less snapshots then observation time steps'

        if return_copy:
            x = self.copy()
        else:
            x = self

        new_list = main_loop(x)

        del x[:]
        x.extend(new_list)

        for A in x:
            self.bool_int_matrix(A)

        if return_copy:
            return x
        else:
            return

    def node_activity_series(self, norm=True):
        """ returns the number of active nodes for each snapshot. """
        active = {}
        n = float(self.number_of_nodes)
        if norm:
            norma = float(n)
        else:
            norma = 1

        for i in range(len(self)):
            outs = sp.coo_matrix(self[i].sum(axis=1))
            ins = sp.coo_matrix(self[i].sum(axis=0))
            nodes1 = set(outs.row)
            nodes2 = set(ins.col)

            nodes = nodes1.union(nodes2)
            active[i] = len(nodes) / norma

        return active

    def edge_activity_series(self, norm=True):
        """ Dict {time:matrix_density} """
        da = {}
        n = float(self.number_of_nodes)
        if norm:
            norma = n * (n-1.0)
        else:
            norma = 1.0

        for i in range(len(self)):
            da[i] = float(self[i].nnz) / norma

        return da

    def shift_start_time(self, new_start_time, return_copy=False):
        """ Returns list of adjacency matrices with new ordering, beginning with
            Index new_start_time using periodic boundary conditions.
        """
        assert new_start_time <= len(self)-1, \
            'new_start_time must be in network observation time.'
        if return_copy:
            x = self.copy()
            x.extend(self[:new_start_time])
            del x[:new_start_time]
            return x
        else:
            self.extend(self[:new_start_time])
            del self[:new_start_time]

    def GST(self, return_copy=False):
        # alias
        self.time_shuffled(return_copy)

    def time_shuffled(self, return_copy=False):
        """ Shuffle times occurence times for each snapshot.

        """
        if return_copy:
            x = self[:]
        else:
            x = self

        random.shuffle(x)

        if return_copy:
            return x
        else:
            return

    def TR(self, return_copy):
        # alias
        self.time_reversed(return_copy)

    def time_reversed(self, return_copy=False):
        """ reverts list and transposes elements

        """
        if return_copy:
            x = self[:]
        else:
            x = self

        for i in range(len(self)):
            x[i] = x[i].transpose()

        x.reverse()

        if return_copy:
            return x
        else:
            return

    def transpose(self, inplace=True):
        """ Transpose all matrices in self.
            If inplace, the object in transposed in place, otherwise a new
            sequence is returned.
        """
        if inplace:
            x = self
        else:
            x = self[:]

        for i in range(len(self)):
            x[i] = x[i].transpose()

        if inplace:
            return
        else:
            return x

    def symmetrize_matrix(self, A):
        """ Returns symmetric version of a non-symm Matrix A as bool-int. """
        M = A + A.transpose()
        M = M.astype('bool')
        M = M.astype('float')
        return M

    def as_undirected(self):
        """ makes every matrix in self symmetric. """
        # if self.is_directed:
        for i in range(len(self)):
            self[i] = self.symmetrize_matrix(self[i])
        self.is_directed = False
        # else:
        #    raise NotImplementedError, "Network is already undirected."

    def clustering_matrix2vector(self, in_file):
        """ Reads file and returns vector from matrix """
        C = mmread(in_file)
        C = lil_matrix(C)
        x = [0.0 for i in range(C.shape[0])]

        indices = zip(C.nonzero()[0], C.nonzero()[1])
        for i, j in indices:
            x[i + j] += C[i, j]

        return x

    def __random_combination(self, iterable, r, with_replacement=False):
        """ Random selection from
            itertools.combinations_with_replacement(iterable, r).

            Parameters
            ----------
            iterable: iterable
                list where samples are drawn from.

            r: int
                number of elements to be sampled

            with_replacement: boolean (optional, default=False)
                if True, combinations with i<=j<=k are returned, if False i<j<k.
        """
        pool = tuple(iterable)
        n = len(pool)
        if with_replacement:
            indices = sorted(random.randrange(n) for i in xrange(r))
        else:
            indices = sorted(random.sample(xrange(n), r))
        return tuple(pool[i] for i in indices)

    def clustering_matrix(self, limit=None, random_iterations=True, replacement=False):
        """ Computes the matrix of clustering coefficients of
            a matrix sequence.

            Parameters
            ----------
            limit: int, optional (default=None)
                Number of time steps to be considered.

            random_iterations: Boolean, optional (default=True)
                If True, sample time triples are considered

            replacement: Boolean, optional (default=False)
                If True, time indices follow the condition i=<j<=k, and i<j<k, if False.

        """
        def triple_product(M1, M2, M3):
            # Product of three matrices
            a3 = M1 * M2 * M3
            tr = (a3.diagonal()).sum()

            clu_norm = (M1 * M2).sum() - ((M1 * M2).diagonal()).sum()
            clu_norm += (M1 * M3).sum() - ((M1 * M3).diagonal()).sum()
            clu_norm += (M2 * M3).sum()-((M2 * M3).diagonal()).sum()

            return tr, clu_norm

        if limit:
            n = limit
        else:
            n = len(self)

        domain = range(n)
        C = lil_matrix((n, n), dtype='float')
        # c=[]

        if random_iterations:
            for l in range(random_iterations):
                (i, j, k) = \
                    self.__random_combination(domain, 3, replacement)
                trace, c_norm = triple_product(self[i], self[j], self[k])
                if c_norm > 0.0:
                    C[j-i, k-j] += float(trace) / c_norm
                    # c.append((i,j,k,float(trace)/c_norm))
        else:
            for (i, j, k) in itertools.combinations(domain, 3):
                trace, c_norm = triple_product(self[i], self[j], self[k])
                if c_norm > 0.0:
                    C[j-i, k-j] += float(trace) / c_norm
                    # c.append((i,j,k,float(trace)/c_norm))

        return C

    def write(self, fname):
        """ writes self to txtfile.
        """
        # generate edge list
        t_edges = []
        for i in range(len(self)):
            # print "extracting edges ",i
            indices = zip(self[i].nonzero()[0], self[i].nonzero()[1])
            to_add = [(u, v, i) for u, v in indices]
            t_edges.extend(to_add)

        # edge list as set for file storage
        t_edges_set = set(t_edges)
        # remove double edges, if undirected
        if not self.is_directed:
            print("removing bidirectional links...")
            for (u, v, d) in t_edges:
                if (v, u, d) in t_edges_set and (u, v, d) in t_edges_set:
                    t_edges_set.remove((v, u, d))

        # write file
        g = file(fname, 'w+')
        for e in t_edges_set:
            wstring = ''
            for j in range(1, len(e)):
                wstring += '\t' + str(e[j])
            g.writelines((str(e[0]) + wstring + '\n'))
        g.close
        return

    def matricesCreation(self):
        """ creates list of sparse matrices from input file """
        edges = loadtxt(self.fname, dtype=int, usecols=self.cols)
        _, _, days = np.array(list(zip(*edges)))

        if not self.first_day:
            self.first_day = min(days)
        if not self.last_day:
            self.last_day = max(days)

        # use only times between firsttime and lasttime
        edges = [(u, v, d) for u, v, d in edges if
                 (d >= self.first_day) and (d <= self.last_day)]

        # get dictionary of new indices and write map-file
        re_dct = self.reindex(edges)
        if self.label_file:
            g = file('oldindex_matrixfriendly.txt', 'w+')
            for k in re_dct:
                g.writelines((str(k) + '\t' + str(re_dct[k]) + '\n'))
            g.close

        # reindex using this dictionary
        edges = [(re_dct[u], re_dct[v], d) for u, v, d in edges]

        edges = self.groupByTime(edges)

        # the actual construction of the sparse matrices
        mx_index = len(re_dct)
        for d, es in edges:
            us = [u for u, v in es]
            vs = [v for u, v in es]
            bs = [True for i in range(len(es))]

            m = csr_matrix((bs, (us, vs)), shape=(mx_index, mx_index),
                           dtype=np.int32)
            self.append(m)

    def bool_int_matrix(self, M):
        """ Returns matrix with only np.int64: ones. """
        M.data = np.ones_like(M.data)

    def unfold_accessibility(self, verbose=True,
                             return_accessibility_matrix=False):
        """ Unfold accessibility storing path density.

        """
        P = self[0].copy()
        D = sp.identity(self.number_of_nodes, dtype=np.int32)
        P = P + D
        cumu = [P.nnz]

        for i in range(1, len(self)):
            if verbose:
                print('unfolding accessibility. Step ', i, 'non-zeros: ', P.nnz)
            self.bool_int_matrix(P)
            try:
                P = P + P * self[i]
            except:
                print('Break at t = ', i)
                break
            cumu.append(P.nnz)
        else:
            print('---> Unfolding complete.')

        if return_accessibility_matrix:
            P = P.astype('bool')
            P = P.astype('int')
            return P, cumu
        else:
            return cumu

    def unfold_accessibility_memory_efficient(self, return_ranges=False):
        """ Computes path density step by step for single nodes.

            Parameters
            ----------
            return ranges: boolean, optional (default=False)
                If True, the method returns a tuple with path density over time
                and the range for every node.

            Returns
            -------
            Returns a numpy vector, where indices are the time steps and values
            are path densities (not normalized).

            If '''return ranges''', returns a tuple with the above path
            denisities and the range of the nodes as a dictionary.

            Usage
            -----
            >>> c = At.unfold_accessibility_memory_efficient()
            >>> c, r = At.unfold_accessibility_memory_efficient(True)

        """
        all_paths = zeros(len(self), dtype=int)
        ranges = {}

        for node in range(self.number_of_nodes):
            print('Computing accessibility for node ', node+1,\
                    ' of ', self.number_of_nodes)
            single_node_SI = self.unfold_accessibility_single_node(node)
            all_paths += single_node_SI
            ranges[node] = single_node_SI[-1]

        if return_ranges:
            return (all_paths, ranges)
        else:
            return all_paths

    def unfold_accessibility_single_node(self, start):
        """ Accessibility of one node. Returns a numpy vector containing
            the number of nonzeros for every timestep.
        """
        # init
        x = sp.coo_matrix(([1], ([0], [start])),
                          shape=(1, self.number_of_nodes), dtype=int)
        x = x.tocsr()

        # these 2 lines are not in the for-loop to be
        # optically consistent with the matrix version.
        x = x + x * self[0]
        cumu = [x.nnz]

        for t in range(1, len(self)):
            x = x + x * self[t]
            cumu.append(x.nnz)

        return np.array(cumu)

    def trace_forward(self, start, stop=None):
        """ same as unfold_accessibility_single_node, but returns all
            nodes reached during traversal.
        """
        if not stop:
            maxtime = len(self)

        # init
        x = sp.coo_matrix(([1], ([0], [start])),
                          shape=(1, self.number_of_nodes), dtype=int)
        x = x.tocsr()

        # these 2 lines are not in the for-loop to be
        # optically consistent with the matrix version.
        x = x + x * self[0]
        cumu = {}

        x = x.tocoo()
        cumu[0] = set(x.col)

        for t in range(1, maxtime):
            x = x + x * self[t]
            x = x.tocoo()
            cumu[t] = set(x.col)

        return cumu

if __name__ == "__main__":
    At = AdjMatrixSequence("../edgelists/sociopatterns_hypertext.dat",
                           directed=True)
    # print len(At), At[0]
    # At.coarse_grain(184)
    # print len(At), At[0]
    # print At.static_path_density(normalize=True)
    # c = At.unfold_accessibility()
