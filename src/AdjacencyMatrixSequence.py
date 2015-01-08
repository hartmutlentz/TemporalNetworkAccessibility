#! /usr/bin/python

from scipy.sparse import coo_matrix, csr_matrix, lil_matrix, dok_matrix
import scipy.sparse as sp
import scipy.stats
import numpy as np
from numpy import loadtxt, zeros, savetxt
import scipy
from random import sample
from collections import defaultdict
from scipy.io import mmread, mmwrite, loadmat, savemat
import random


class AdjMatrixSequence(list):
    """
    list of sparse matrixes

    class inherits from list.
    the constructor expects filename of u,v,w,d-edgelist

    it constructs a list, where each entry is
    a sparse matrix (default: csr)
    the matrices in the list are ordered by day

    the indices of the matrix represent the nodes, but
    nodes are reindexed to [0..number_of_nodes]
    """
    def __init__(self, edgelist_fname, directed, write_label_file=True,
                 columns=(0, 1, 2), firstday=None, lastday=None):
        list.__init__(self)
        #if edgelist_fname == fs.dataPath("D_sf_uvwd_cmpl.txt"):
        #    self.first_day = 2555 #2008-2010
        #else:
        #    self.first_day = 0#1825 #2006-2010
        #self.last_day = 9#3650
        self.first_day = firstday
        self.last_day = lastday
        self.fname = edgelist_fname
        self.cols = columns
        self.label_file = write_label_file
        self.is_directed = directed

        self.matricesCreation()
        if not self.is_directed:
            self.as_undirected()
        self.number_of_nodes = scipy.shape(self[0])[0]

    def __scipy_version_for_large_matrices(self, ref="0.14.0"):
        """ Checks if the installed version of Scipy is at least ref.
            This is necessary for handling very large sparse matrices.
            For scipy versions < 0.14.0. the variables indptr and indices are
            stored as int32. Thus te number of non-zero entries is restricted to
            2^31 \sim 10^9 elements.
            Returns True, if installed Version > ref. False otherwise.
        """
        # split versions into numbers
        x = scipy.__version__.split(".")
        y = ref.split(".")
        
        # compare adjusted lengths
        if len(x)==len(y): return x >= y
        elif len(x) < len(y): return x >= y[:len(x)]
        else: return x[:len(y)] >= y
    
    def info_scipy_version(self):
        """ Print information about scipy version and maximum Matrix size """
        if self.__scipy_version_for_large_matrices():
            print ("Scipy version can handle matrices with more than " +
                   "2^31 nonzero elements.")
        else:
            print "Number of nonzero elements is restricted to 2^31."

    def groupByDays(self, edges):
        """ returns list of tupels: [(d,[(u,v),...]),...] """
        dct = defaultdict(list)
        for u, v, d in edges:
            dct[d].append((u, v))
        dct_s = dict.fromkeys(range(0, self.last_day-self.first_day), [])
        for d in dct:
            dct_s[d-self.first_day] = dct[d]
        return dct_s.items()

    def reindex(self, edges):
        """ for any index in edges returns dict with new_indices
            [0...number_of_unique_indices]
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
            print twindow
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
            return float(nlinks) / (float((n - 2) * (n**2 - n))
                                    + float(n**2 - n))
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

    def LCCs(self):
        """ returns information about the size of the LCC and the mean degree
            of the snapshots.
        """
        crit = {}
        for i, M in enumerate(self):
            print "LCC ", i
            crit[i] = (self.__matrix_mean_degree(M), self.__matrix_LCC_size(M))

        return crit

    def path_density_of_A(self, A):
        """ The path density of an Adjacency Matrix A """
        paths = 0
        n = scipy.shape(A)[0]

        for i in range(n):
            out_size =\
                len(sp.csgraph.depth_first_order(A,
                    i, return_predecessors=False))-1
            paths += out_size

        return float(paths) / n**2

    def static_path_density(self):
        """ Returns dict. {Aggregation depth: static path density}

        """
        pd = []
        for i, Cn in enumerate(self.step_by_step_aggregation()):
            print 'Static path density. Step ', i
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

    def cumulated(self, start=0, ende=None):
        """ Returns Cumulated Graph as Matrix """
        C = csr_matrix((self.number_of_nodes, self.number_of_nodes),
                       dtype=np.int32)

        if ende:
            e = ende
        else:
            e = len(self)

        for i in range(start, e):
            C = C + self[i]

        return C

    def daily_activity(self):
        """ Dict {day:matrix_density} """
        da = {}
        n = float(self.number_of_nodes)
        norma = n * (n-1.0)

        for i in range(len(self)):
            da[i] = float(self[i].nnz) / norma

        return da

    def GST(self, return_copy=False):
        # alias
        self.time_shuffled(return_copy)

    def time_shuffled(self, return_copy=False):
        """

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

    def symmetrize_matrix(self, A):
        """ Returns symmetric version of a non-symm Matrix A as bool-int. """
        M = A + A.transpose()
        M = M.astype('bool')
        M = M.astype('float')
        return M

    def as_undirected(self):
        """ makes every matrix in self symmetric. """
        #if self.is_directed:
        for i in range(len(self)):
            self[i] = self.symmetrize_matrix(self[i])
        self.is_directed = False
        #else:
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

    def __random_combination_with_replacement(self, iterable, r):
        """ Random selection from
            itertools.combinations_with_replacement(iterable, r)
        """
        pool = tuple(iterable)
        n = len(pool)
        indices = sorted(random.randrange(n) for i in xrange(r))
        return tuple(pool[i] for i in indices)

    def clustering_matrix(self, limit=None, random_iterations=False):
        """ Computes the matrix of clustering coefficients of
            a matrix sequence.
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
        #c=[]

        if random_iterations:
            for l in range(random_iterations):
                (i, j, k) = \
                    self.__random_combination_with_replacement(domain, 3)
                trace, c_norm = triple_product(self[i], self[j], self[k])
                if c_norm > 0.0:
                    C[j-i, k-j] += float(trace) / c_norm
                    #c.append((i,j,k,float(trace)/c_norm))
        else:
            for (i, j, k) in itertools.combinations(domain, 3):
                trace, c_norm = triple_product(self[i], self[j], self[k])
                if c_norm > 0.0:
                    C[j-i, k-j] += float(trace) / c_norm
                    #c.append((i,j,k,float(trace)/c_norm))

        return C
        #return c

    def write(self, fname):
        """ writes self to txtfile.
            If network is undirected, edge-pairs appear twice.
        """
        # generate edge list
        t_edges = []
        for i in range(len(self)):
            #print "extracting edges ",i
            indices = zip(self[i].nonzero()[0], self[i].nonzero()[1])
            to_add = [(u, v, i) for u, v in indices]
            t_edges.extend(to_add)

        # edge list as set for file storage
        t_edges_set = set(t_edges)
        # remove double edges, if undirected
        if not self.is_directed:
            print "removing bidirectional links..."
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

        # first and last days
        _, _, days = loadtxt(self.fname, dtype=int, usecols=self.cols,
                             unpack=True)
        if not self.first_day:
            self.first_day = min(days)
        if not self.last_day:
            self.last_day = max(days)

        # use only days between FIRSTDAY and LASTDAY
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

        edges = self.groupByDays(edges)

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
        M = M.astype('bool')
        M = M.astype('i')

    def unfold_accessibility(self, return_accessibility_matrix=False):
        """ Unfold accessibility storing path density.

        """
        P = self[0].copy()
        D = sp.identity(self.number_of_nodes, dtype=np.int32)
        P = P + D
        cumu = [P.nnz]

        for i in range(1, len(self)):
            print 'unfolding accessibility. Step ', i, 'non-zeros: ', P.nnz
            self.bool_int_matrix(P)
            cumu.append(P.nnz)
            try:
                P = P + P * self[i]
            except:
                print 'Break at t = ', i
                break
        else:
            print '---> Unfolding complete.'

        if return_accessibility_matrix:
            P = P.astype('bool')
            P = P.astype('int')
            return P, cumu
        else:
            return cumu


if __name__ == "__main__":
    from pprint import pprint

    At = AdjMatrixSequence("sociopatterns_hypertext.dat", directed=False)
    print len(At)
    c = At.unfold_accessibility()
