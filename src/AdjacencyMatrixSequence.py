#! /usr/bin/python
"""

Provides Class AdjMatrixSequence for temporal network analysis.

Networks are represented as a sequence of adjacency matrices and each matrix
is a snapshot of the network.

Author: Hartmut H. K. Lentz
https://github.com/hartmutlentz/TemporalNetworkAccessibility
"""

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
import sys
from collections import deque


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

        self.create_matrices()
        if not self.is_directed:
            self.as_undirected()
        self.number_of_nodes = np.shape(self[0])[0]
        self.check_py_version()

    @staticmethod
    def check_py_version():
        """Check version."""
        assert sys.version_info > (3,), ("You are using python 2. Please use "
                                         "python 3.\nPython 3.6 or greater is "
                                         "recommended.")

    def copy(self):
        """Alias."""
        return copy.copy(self)

    def deepcopy(self):
        """Alias."""
        return copy.deepcopy(self)

    def __str__(self):
        """
        Magic method for printing.

        Returns
        -------
        str
            description.

        """
        a = "Temporal network stored as AdjMatrixSequence.\n"
        b = "Number of nodes: " + str(self.number_of_nodes) + "; "
        c = "Number of edges: " + str(self.get_number_of_edges()) + "; "
        d = "Time steps: " + str(self.last_day - self.first_day)

        return a + b + c + d

    @staticmethod
    def __scipy_version_for_large_matrices(ref="0.14.0"):
        """
        Check if the installed version of Scipy is at least ref.

        This is necessary for handling very large sparse matrices.
        For scipy versions < 0.14.0. the variables indptr and indices are
        stored as int32. Thus, the number of non-zero entries is restricted to
        2^31 ~ 10^9 elements.
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
        """Print information about scipy version and maximum Matrix size."""
        if self.__scipy_version_for_large_matrices():
            print("Scipy version can handle matrices with more than " +
                  "2^31 nonzero elements.")
        else:
            print("Number of nonzero elements is restricted to 2^31.")

    def groupbytime(self, edges):
        """Return list of tupels: [(d,[(u,v),...]),...]."""
        dct = defaultdict(list)
        for u, v, d in edges:
            dct[d].append((u, v))
        dct_s = dict.fromkeys(range(0, self.last_day-self.first_day), [])
        for d in dct:
            dct_s[d-self.first_day] = dct[d]
        return dct_s.items()

    def reindex(self, edges):
        """
        For any index in edges returns dict with new_indices.

        [0...number_of_unique_indices].
        """
        us, vs, ds = zip(*edges)
        nodes = set(us) | set(vs)
        old_to_new = {}
        for i, n in enumerate(nodes):
            old_to_new[n] = i
        return old_to_new

    def all_time_windows(self, max_window=None):
        """Summation over all time windows."""
        p_corr = {}
        if not max_window:
            mw = len(self)
        else:
            mw = max_window

        for twindow in range(mw):
            print("All time windows:", twindow)
            p_corr[twindow] = self.single_time_window(twindow)
        return p_corr

    def single_time_window(self, windowsize):
        """Summation over a time window."""
        prop_corr = []

        for i in range(len(self) - windowsize):
            prop_corr.append(self.two_link_density(i, i+windowsize))

        return (scipy.stats.scoreatpercentile(prop_corr, 25),
                scipy.stats.scoreatpercentile(prop_corr, 50),
                scipy.stats.scoreatpercentile(prop_corr, 75))

    def get_number_of_edges(self):
        """
        Return total number of edges.

        Returns
        -------
        int
            number of edges.

        """
        edges = [A.nnz for A in self]

        return sum(edges)

    def two_link_density(self, index1, index2, norm=True):
        """Link density for 2 step paths."""
        C = self[index1] * self[index2]

        nlinks = C.sum()
        n = self.number_of_nodes

        if norm:
            return float(nlinks) / (float((n - 2) * (n**2 - n)) +
                                    float(n**2 - n))
        else:
            return float(nlinks) / float((n**2 - n))

    def deep_product(self, twindow=1, start=0):
        """Product A_1*A_7*A_14... ."""
        C = self[start].copy()
        links = {start: (C * C).sum()}

        for i in range(start+twindow, len(self)-twindow, twindow):
            C = C * self[i]
            links[i] = C.nnz
            if C.nnz == 0:
                break

        return links

    def __matrix_mean_degree(self, A):
        """Mean degree of A."""
        return float(A.nnz) / (self.number_of_nodes - 1)

    def __matrix_LCC_size(self, A):
        """
        Size of the L(S)CC of a matrix.

        Note that the matrices here are interpreted as (bi)directed nets.
        """
        n, m = sp.csgraph.connected_components(A, connection='strong')
        lcc_size = np.max(np.bincount(m))
        if lcc_size == 1:
            return 0.0
        else:
            return lcc_size / float(self.number_of_nodes)

    @staticmethod
    def average_path_length(M, diameter):
        """
        Return average pathlength of a snapshot.

        All path lengths > 1 are considered.
        """
        x = [(M**i).nnz for i in range(2, diameter+1)]

        return np.mean(x)

    def long_paths_per_snapshot(self, max_path_length):
        """
        Compute the number of paths longer than 1 in each snapshot.

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
        """
        Replace A_i by sum _{i=1} ^diameter A_i.

        This takes into account paths of length > 1 in each snapshot.
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

        print("---> paths up to length ", diameter,
              " are now considered in snapshots.")

        return

    def LCCs(self):
        """
        Return information about the size of the LCC.

        ... and the mean degree of the snapshots.
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
    @staticmethod
    def LSCC_nodes(A):
        """
        Return nodes of the largest stringly connected component.

        Network is given by a sparse matrix A.

        Parameters
        ----------
        A - scipy.sparse matrix
            Adjacency matrix of the network. Assumed to be directed.

        Returns
        -------
        set
            Node labels of the nodes in the LSCC.
        """
        n, labs = sp.csgraph.connected_components(A, connection="strong")

        nodes = set()
        for i in range(len(labs)):
            if labs[i] == 0:
                nodes.add(i)

        return nodes

    def path_density_of_A(self, A, normalize=False):
        """
        Return the path density of a network.

        Network is given by adjacency matrix A.
        Selfloops are explicitely included.
        """
        n = A.shape[0]  # scipy.shape(A)[0]
        all_nodes = set(range(n))

        lscc_nodes = self.LSCC_nodes(A)
        rest_nodes = all_nodes - lscc_nodes

        paths = 0

        # LSCC nodes
        lscc_sample_node = next(iter(lscc_nodes))
        lscc_range = len(sp.csgraph.depth_first_order(A,
                         lscc_sample_node, return_predecessors=False))
        for _ in lscc_nodes:
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
        Path density of the aggregated network.

        This algorithm first
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
        >>> # This method works in an At-Object only.
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

    def step_by_step_static_path_density(self, ende=None, verbose=False):
        """
        Return list.

        [index=Aggregation depth: static path density]

        """
        verboseprint = print if verbose else lambda *a, **k: None

        pd = []
        for i, Cn in enumerate(self.step_by_step_aggregation(ende)):
            verboseprint('Static path density. Step ', i)
            pd.append(self.path_density_of_A(Cn))

        return pd

    def step_by_step_aggregation(self, ende=None):
        """
        Return matrix list of all aggregated networks.

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
        """Return Cumulated Graph as Matrix."""
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
        """
        Coarse grain the list, i.e. partial aggregation of the network.

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
        """Return the number of active nodes for each snapshot."""
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
        """Return dict {time:matrix_density}."""
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
        """
        Return list of adjacency matrices with new ordering.

        beginning with Index new_start_time using periodic boundary conditions.
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

    @staticmethod
    def random_submatrix(A, p=0.5):
        """
        Return a random subset of a sparse matrix.

        Dimension is not changed.
        Output-values are Boolean, i.e. 1 (int).

        Parameters
        ----------
        A : scipy sparse matrix
            An adjacency matrix or any other sparse matrix.
        p : float, optional
            probability that an entry remains in the matrix. Thus, 1-p is the
            removal probability. The default is 0.5.

        Returns
        -------
        B : scipy sparse csr matrix
            A matrix containing a random subset of the input. Dimension is not
            changed.

        """
        assert 0.0 <= p <= 1.0, "p is a probability and must be\
            0 <= p <= 1."

        indices = np.column_stack(A.nonzero())
        number_of_samples = round(indices.shape[0] * p)

        index_sample = np.random.choice(indices.shape[0], number_of_samples,
                                        replace=False)
        rows_and_columns = indices[index_sample, :]

        row, col = rows_and_columns.T
        data = np.ones(len(row), dtype=int)
        B = sp.csr_matrix((data, (row, col)), shape=A.shape)

        return B

    def dilute(self, p=0.5):
        """
        Remove edges from the network randomly with probability 1-p.

        Thus, p is the probablity that an edge remains in the network.
        CAUTION: In-place operation. Make a copy first.

        Parameters
        ----------
        p : float, optional
            Probability for an edge to remain in the network.
            The default is 0.5.

        Returns
        -------
        None.

        """
        for i in range(len(self)):
            self[i] = self.random_submatrix(self[i], p)

    def GST(self, return_copy=False):
        """Alias."""
        if return_copy:
            return self.time_shuffled(return_copy)
        else:
            self.time_shuffled(return_copy)

    def time_shuffled(self, return_copy=False):
        """Shuffle times occurence times for each snapshot."""
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
        """Alias."""
        if return_copy:
            return self.time_reversed(return_copy)
        else:
            self.time_reversed(return_copy)

    def time_reversed(self, return_copy=False):
        """Revert list and transposes elements."""
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
        """
        Transpose all matrices in self.

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

    @staticmethod
    def symmetrize_matrix(A):
        """Return symmetric version of a non-symm Matrix A as bool-int."""
        M = A + A.transpose()
        M = M.astype('bool')
        M = M.astype('float')
        return M

    def as_undirected(self):
        """Make every matrix in self symmetric."""
        # if self.is_directed:
        for i in range(len(self)):
            self[i] = self.symmetrize_matrix(self[i])
        self.is_directed = False
        # else:
        #    raise NotImplementedError, "Network is already undirected."

    @staticmethod
    def clustering_matrix2vector(in_file):
        """Read file and returns vector from matrix."""
        C = mmread(in_file)
        C = lil_matrix(C)
        x = [0.0 for _ in range(C.shape[0])]

        indices = zip(C.nonzero()[0], C.nonzero()[1])
        for i, j in indices:
            x[i + j] += C[i, j]

        return x

    @staticmethod
    def __random_combination(iterable, r, with_replacement=False):
        """
        Random selection.

        from itertools.combinations_with_replacement(iterable, r).

        Parameters
        ----------
        iterable: iterable
            list where samples are drawn from.

        r: int
            number of elements to be sampled

        with_replacement: boolean (optional, default=False)
            if True, combinations with i<=j<=k are returned,
            if False i<j<k.
        """
        pool = tuple(iterable)
        n = len(pool)
        if with_replacement:
            indices = sorted(random.randrange(n) for _ in range(r))
        else:
            indices = sorted(random.sample(range(n), r))
        return tuple(pool[i] for i in indices)

    def clustering_matrix(self, limit=None, random_iterations=True,
                          replacement=False):
        """
        Compute the matrix of clustering coefficients of a matrix sequence.

        Parameters
        ----------
        limit: int, optional (default=None)
            Number of time steps to be considered.

        random_iterations: Boolean, optional (default=True)
            If True, sample time triples are considered

        replacement: Boolean, optional (default=False)
            If True, time indices follow the condition i <= j <= k, and
            i < j < k, if False.

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
            for _ in range(random_iterations):
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
        """Write self to txtfile."""
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
        g = open(fname, 'w+')
        for e in t_edges_set:
            wstring = ''
            for j in range(1, len(e)):
                wstring += '\t' + str(e[j])
            g.writelines((str(e[0]) + wstring + '\n'))
        g.close()
        return

    def create_matrices(self):
        """Create list of sparse matrices from input file."""
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
            g = open('oldindex_matrixfriendly.txt', 'w+')
            for k in re_dct:
                g.writelines((str(k) + '\t' + str(re_dct[k]) + '\n'))
            g.close()

        # reindex using this dictionary
        edges = [(re_dct[u], re_dct[v], d) for u, v, d in edges]

        edges = self.groupbytime(edges)

        # the actual construction of the sparse matrices
        mx_index = len(re_dct)
        for d, es in edges:
            us = [u for u, v in es]
            vs = [v for u, v in es]
            bs = [True for _ in range(len(es))]

            m = csr_matrix((bs, (us, vs)), shape=(mx_index, mx_index),
                           dtype=np.int32)
            self.append(m)

    @staticmethod
    def bool_int_matrix(M):
        """Return matrix with only np.int64: ones."""
        M.data = np.ones_like(M.data)

    def si_model(self, p=0.1, verbose=True, return_accessibility_matrix=False):
        """
        Compute an SI-model (susceptible, infected) on the network.

        This is statistically equivalent to dilute the network first,
        and then unfold the accessibility matrix.

        Parameters
        ----------
        p : float, optional
            Infection probability per temporal edge. Default is 0.1.

        verbose : Boolean, optional
            For text output. The default is True.

        return_accessibility_matrix : Boolean, optional
            Returns the whole accessibility matrix. The matrix can be huge for
            large networks. The default is False.

        Returns
        -------
        list
            Cumulative path density vs. time.

        """
        assert p <= 1.0, "Probability must be <= 1."
        assert p >= 0.0, "Probability must be >= 0."
        verboseprint = print if verbose else lambda *a, **k: None

        P = self[0].copy()
        P.multiply(csr_matrix((
            np.random.random_sample(P.data.shape) < p, P.indices, P.indptr),
            shape=P.shape))

        D = sp.identity(self.number_of_nodes, dtype=np.int32)
        P = P + D
        cumu = [P.nnz]

        for i in range(1, len(self)):
            verboseprint(i, end=" ")
            self.bool_int_matrix(P)
            try:
                X = self[i].multiply(csr_matrix((
                        np.random.random_sample(self[i].data.shape) < p,
                        self[i].indices, self[i].indptr),
                        shape=P.shape))
                P = P + P * X
            except (KeyboardInterrupt, SystemExit, MemoryError):
                print('\nBreak at t = ', i)
                break
            cumu.append(P.nnz)
        else:
            print('\n---> Unfolding complete.')

        if return_accessibility_matrix:
            P = P.astype('bool')
            P = P.astype('int')
            return P, cumu
        else:
            return cumu

    def unfold_accessibility(self, verbose=True,
                             return_accessibility_matrix=False):
        """
        Compute path density for all nodes.

        Parameters
        ----------
        verbose : Boolean, optional
            For text output. The default is True.
        return_accessibility_matrix : Boolean, optional
            Returns the whole accessibility matrix. The matrix can be huge for
            large networks. The default is False.

        Returns
        -------
        list
            Cumulative path density vs. time.

        Usage
        -----
        >>> c = At.unfold_accessibility()
        >>> c, r = At.unfold_accessibility_memory_efficient(
                    return_accessibility_matrix=True)

        """
        verboseprint = print if verbose else lambda *a, **k: None

        P = self[0].copy()
        D = sp.identity(self.number_of_nodes, dtype=np.int32)
        P = P + D
        cumu = [P.nnz]

        for i in range(1, len(self)):
            verboseprint(i, end=" ")
            self.bool_int_matrix(P)
            try:
                P = P + P * self[i]
            except (KeyboardInterrupt, SystemExit, MemoryError):
                print('\nBreak at t = ', i)
                break
            cumu.append(P.nnz)
        else:
            print('\n---> Unfolding complete.')

        if return_accessibility_matrix:
            P = P.astype('bool')
            P = P.astype('int')
            return P, cumu
        else:
            return cumu

    def unfold_accessibility_memory_efficient(self, return_ranges=False):
        """
        Compute path density step by step for single nodes.

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
            print('Computing accessibility for node ', node+1,
                  ' of ', self.number_of_nodes)
            single_node_SI = self.unfold_accessibility_single_node(node)
            all_paths += single_node_SI
            ranges[node] = single_node_SI[-1]

        if return_ranges:
            return all_paths, ranges
        else:
            return all_paths

    def unfold_accessibility_single_node(self, start):
        """
        Compute accessibility of one node.

        Returns a numpy array containing the number of nonzeros
        for every timestep.
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

    def unfold_accessibility_multi_nodes(self, start):
        """
        Compute accessibility of multiple, but not all, nodes.

        Returns a numpy array containing the number of nonzeros
        for every timestep.
        """
        # init
        row = np.zeros(len(start))
        col = np.array(start)
        data = np.ones(len(start))

        x = sp.coo_matrix((data, (row, col)),
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

    def unfold_accessibility_with_sentinels(self, sentinels,
                                            start_node=None,
                                            start_time=0,
                                            stop_at_detection=False):
        """
        Unfold the accessibility graph including sentinel nodes.

        Sentinel nodes are for disease detection. The method returns the
        arrival times of an epidemic at the *sentinels* starting at
        *start_node*. In order to simulate SI-like spreading the network
        should be diluted.

        Parameters
        ----------
        sentinels : list
            list of sentinel nodes.
        start_node : int, optional
            index node of the epidemic. The default is None. If default is used
            staring node is chosen at random.
        start_time : int, optional
            starting time for the epidemic. The default is 0.
        stop_at_detection : Boolean, optional
            If true, the epidemic is stopped, when it arrives at any sentonel
            node. The default is False.

        Returns
        -------
        Dictionary. Keys are sentinel nodes and values are tuples with
        (arrival time, outbreak size).

        """
        # raise error if sentinel node not in network
        if max(sentinels) >= self.number_of_nodes:
            raise ValueError("Sentinel node not in network.")

        # set start_node for epidemic
        if start_node or start_node == 0:
            start = start_node
        else:
            start = np.random.randint(self.number_of_nodes)
        print("Starting epidemic at node ", start)

        # state array
        x = sp.csr_matrix(([1], ([0], [start])),
                          shape=(1, self.number_of_nodes), dtype=int)

        # sentinels array
        row = np.zeros(len(sentinels))
        col = np.array(sentinels)
        data = np.ones(len(sentinels))
        sen_nodes = sp.csr_matrix((data, (row, col)),
                                  shape=(1, self.number_of_nodes), dtype=int)

        # sentinel arrival time dict
        arrival_times = dict()

        if stop_at_detection:
            for t in range(start_time, len(self)):
                x = x + x * self[t]
                if (x.multiply(sen_nodes)).nnz > 0:
                    infected_sentinels = set((x.multiply(sen_nodes) != 0)
                                             .nonzero()[1])
                    arrival_times.update({node: (t, x.nnz) for node in
                                          infected_sentinels})
                    break
        else:
            for t in range(start_time, len(self)):
                x = x + x * self[t]
                infected_sentinels = set((x.multiply(sen_nodes) != 0)
                                         .nonzero()[1])
                new_infected = infected_sentinels - arrival_times.keys()
                arrival_times.update({node: (t, x.nnz) for node in
                                      new_infected})

        return arrival_times

    def unfold_accessibility_sir_constant_recovery(self, p_si=0.5, recovery=10, return_accessibility_matrix = False):

        """
        Compute path density for all nodes. Infected nodes will recover after a given timespan.
        Recovered nodes can not be infected again.
        Parameters
        ----------
        p_si : float, optional
            probability of becoming infected from contact with an infected node. Default is 0.5.
        recovery : int
            Time after which an infected node recovers and cannot be infected again.
        return_accessibility_matrix : Boolean, optional
            Returns the whole accessibility matrix. The matrix can be huge for
            large networks. The default is False.
        Returns
        -------
        cumu : list
        a list containing the number of infected nodes for every timestep.
        cumu_rec : list
        a list containing the number of removed nodes for every timestep.
        
        """

        P = self[0].copy()
        D = sp.identity(self.number_of_nodes, dtype=np.int32)
        # the si part
        P = D + D * csr_matrix((np.random.random_sample(self[0].data.shape)<p_si, self[0].indices, self[0].indptr), shape=self[0].shape)
        # the ir part
        empty = csr_matrix((self.number_of_nodes, self.number_of_nodes), dtype=np.int32)
        history = deque([empty for i in range(recovery)], maxlen=recovery)
        history.append(P)
        R = empty
        # save state
        cumu = [P.nnz]
        cumu_rec = [0]

        for t in range(1, len(self)): 
            # the si part
            P = P + P * csr_matrix((np.random.random_sample(self[t].data.shape)<p_si, self[t].indices, self[t].indptr), shape=self[t].shape)
            # the ir part
            R += history[0]
            P -= P.multiply(R.astype("bool"))
            history.append(P)
            # save state
            cumu.append(P.nnz)
            cumu_rec.append(R.nnz)


        if return_accessibility_matrix:
            P += R
            self.bool_int_matrix(P)
            P = P.astype('bool')
            P = P.astype('int')
            return P, cumu, cumu_rec
        else:
            return cumu, cumu_rec

    def unfold_accessibility_sir_random_recovery(self, p_si=0.5, p_ir=0.01, return_accessibility_matrix = False):

        """
        Compute path density for all nodes. Infected nodes will recover spontaneously with a given probability.
        Recovered nodes can not be infected again.
        Parameters
        ----------
        p_si : float, optional
            probability of becoming infected from contact with an infected node. Default is 0.5.
        p_ir : float
            Probability that a node recovers spontaneously. Should be between 0 and 1. Default is 0.01.
        return_accessibility_matrix : Boolean, optional
            Returns the whole accessibility matrix. The matrix can be huge for
            large networks. The default is False.
        Returns
        -------
        cumu : list
        a list containing the number of infected nodes for every timestep.
        cumu_rec : list
        a list containing the number of removed nodes for every timestep.
        
        """

        P = self[0].copy()
        D = sp.identity(self.number_of_nodes, dtype=np.int32)
        # the si part
        P = D + D * csr_matrix((np.random.random_sample(self[0].data.shape)<p_si, self[0].indices, self[0].indptr), shape=self[0].shape)
        #the ir part
        R = P.multiply(csr_matrix((np.random.random_sample(P.data.shape)<p_ir, P.indices, P.indptr), shape=P.shape))
        P -= P.multiply(R.astype("bool"))
        # save state
        cumu = [P.nnz]
        cumu_rec = [R.nnz]

        for t in range(1, len(self)):
            # the si part
            P = P + P * csr_matrix((np.random.random_sample(self[t].data.shape)<p_si, self[t].indices, self[t].indptr), shape=self[t].shape)
            # the ir part
            R += P.multiply(csr_matrix((np.random.random_sample(P.data.shape)<p_ir, P.indices, P.indptr), shape=P.shape))
            P -= P.multiply(R.astype("bool"))
            # save state
            cumu.append(P.nnz)
            cumu_rec.append(R.nnz)

        if return_accessibility_matrix:
            P += R
            P = P.astype('bool')
            P = P.astype('int')
            return P, cumu, cumu_rec
        else:
            return cumu, cumu_rec

    def unfold_accessibility_sir_constant_recovery_single_node(self, p_si=0.5, recovery=10, start_node = None, return_accessibility_matrix = False):
        
        """
        Compute path density for a single node. Infected nodes will recover after a given timespan.
        Recovered nodes can not be infected again.
        Parameters
        ----------
        p_si : float, optional
            probability of becoming infected from contact with an infected node. Default is 0.5.
        recovery : int
            Time after which an infected node recovers and cannot be infected again.
        start_node : int, optional
            First infected node. If None, a random node is chosen.
        return_accessibility_matrix : Boolean, optional
            Returns the whole accessibility matrix. The matrix can be huge for
            large networks. The default is False.
        Returns
        -------
        cumu : list
        a list containing the number of infected nodes for every timestep.
        cumu_rec : list
        a list containing the number of removed nodes for every timestep.
        
        """

        # set start_node for epidemic
        if start_node or start_node == 0:
            start = start_node
        else:
            start = np.random.randint(self.number_of_nodes)
        print("Starting epidemic at node ", start)

        x = sp.csr_matrix(([1], ([0], [start])), shape=(1, self.number_of_nodes), dtype=int)
        # the si part
        x = x + x * csr_matrix((np.random.random_sample(self[0].data.shape)<p_si, self[0].indices, self[0].indptr), shape=self[0].shape)
        # the ir part
        empty = csr_matrix((1, self.number_of_nodes), dtype=np.int32)
        history = deque([empty.copy() for i in range(recovery)], maxlen=recovery)
        R = empty.copy()
        #save state
        cumu = [x.nnz]
        cumu_rec = [0]
        
        for t in range(1, len(self)): 
            # the si part
            x = x + x * csr_matrix((np.random.random_sample(self[t].data.shape)<p_si, self[t].indices, self[t].indptr), shape=self[t].shape)
            # the ir part
            R += history[0]
            x -= x.multiply(R.astype("bool"))
            history.append(x)
            # save state
            cumu.append(x.nnz)
            cumu_rec.append(R.nnz)

        return cumu, cumu_rec
    
    def unfold_accessibility_sir_random_recovery_single_node(self, p_si=0.5, p_ir=0.01, start_node = None, return_accessibility_matrix = False):
        
        """
        Compute path density for a single node. Infected nodes will recover spontaneously with a given probability.
        Recovered nodes can not be infected again.
        Parameters
        ----------
        p_si : float, optional
            probability of becoming infected from contact with an infected node. Default is 0.5.
        p_ir : float
            Probability that a node recovers spontaneously. Should be between 0 and 1. Default is 0.01.
        start_node : int, optional
            First infected node. If None, a random node is chosen.
        return_accessibility_matrix : Boolean, optional
            Returns the whole accessibility matrix. The matrix can be huge for
            large networks. The default is False.
        Returns
        -------
        cumu : list
        a list containing the number of infected nodes for every timestep.
        cumu_rec : list
        a list containing the number of removed nodes for every timestep.
        
        """

        # set start_node for epidemic
        if start_node or start_node == 0:
            start = start_node
        else:
            start = np.random.randint(self.number_of_nodes)
        print("Starting epidemic at node ", start)

        x = sp.csr_matrix(([1], ([0], [start])), shape=(1, self.number_of_nodes), dtype=int)
        # the si part
        x = x + x * csr_matrix((np.random.random_sample(self[0].data.shape)<p_si, self[0].indices, self[0].indptr), shape=self[0].shape)
        # the ir part
        R = x.multiply(csr_matrix((np.random.random_sample(x.data.shape)<p_ir, x.indices, x.indptr), shape=x.shape))
        x -= x.multiply(R.astype("bool"))
        # save state
        cumu = [x.nnz]
        cumu_rec = [R.nnz]

        
        for t in range(1, len(self)): 
            # the si part
            x = x + x * csr_matrix((np.random.random_sample(self[t].data.shape)<p_si, self[t].indices, self[t].indptr), shape=self[t].shape)
            # the ir part
            R += x.multiply(csr_matrix((np.random.random_sample(x.data.shape)<p_ir, x.indices, x.indptr), shape=x.shape))
            x -= x.multiply(R.astype("bool"))
            # save state
            cumu.append(x.nnz)
            cumu_rec.append(R.nnz)

        return cumu, cumu_rec

    def trace_forward(self, start_node, start_time=0, stop=None):
        """
        Similar to unfold_accessibility_single_node.

        But returns all nodes reached during traversal.
        """
        if not stop:
            stop = len(self)

        # init
        x = sp.coo_matrix(([1], ([0], [start_node])),
                          shape=(1, self.number_of_nodes), dtype=int)
        x = x.tocsr()

        # these 2 lines are not in the for-loop to be
        # optically consistent with the matrix version.
        # x = x + x * self[0]
        cumu = {}

        # x = x.tocoo()
        # cumu[0] = set(x.col)

        for t in range(start_time, stop):
            x = x + x * self[t]
            x = x.tocoo()
            cumu[t] = set(x.col)

        return cumu

    def trace_forward_multiple_sources(self, start, stop=None):
        """
        Tracing forward for multiple starting nodes.

        Parameters
        ----------
        start : iterable

        stop (optional) : stop time, default = infty
        """
        if not stop:
            maxtime = len(self)
        else:
            maxtime = stop

        # init
        row = np.zeros(len(start))
        col = np.array(start)
        data = np.ones(len(start))

        x = sp.coo_matrix((data, (row, col)),
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
    print("===== Testing Module AdjacencyMatrixSequence =====\n")
    # the_file = '../edgelists/Test.dat'
    the_file = "../edgelists/sociopatterns_hypertext.dat"
    At = AdjMatrixSequence(the_file, directed=False, write_label_file=False)

    # compute accessibility
    c = At.unfold_accessibility(return_accessibility_matrix=False)
    h = np.gradient(c)

    # Causal fidelity
    causal_paths = c[-1]
    static_paths = At.static_path_density()
    print("---> Causal fidelity is ", float(causal_paths)/float(static_paths))

    # full causal fidelity
    c2 = At.step_by_step_static_path_density()
    c_ff = [c[i]/c2[i] for i in range(len(c))]
    print("---> Functional form of causal fidelity computed")

    At.info_scipy_version()

    # x = At.all_time_windows()
    x = At.deep_product()
    x = At.long_paths_per_snapshot(3)
    At.long_path_correction(3)
    x = At.LCCs()
    x = At.static_path_density()
    x = At.step_by_step_static_path_density()
    x = At.step_by_step_aggregation()
    x = At.cumulated()
    x = At.coarse_grain(2)
    x = At.node_activity_series()
    x = At.edge_activity_series()
    x = At.shift_start_time(0)
    x = At.GST()
    x = At.time_reversed(True)
    x = At.transpose()
    At.as_undirected()
    x = At.clustering_matrix(replacement=True)
    x = At.unfold_accessibility_memory_efficient()
    x = At.unfold_accessibility_single_node(0)
    x = At.trace_forward(0)

    x = At.GST(return_copy=False)

    print(At)
    print(x)
    print("===== Test for AdjacencyMatrixSequence successful. =====")
    # print len(At), At[0]
    # At.coarse_grain(184)
    # print len(At), At[0]
    # print At.static_path_density(normalize=True)
    # c = At.unfold_accessibility()
