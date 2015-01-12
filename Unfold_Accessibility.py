#!/usr/bin/env python
#
#
import sys
sys.path.append('./src')
from AdjacencyMatrixSequence import AdjMatrixSequence
from TemporalNetworkEdgeList import TemporalEdgeList
import Tools


# import an edgelist as sequence of adjacency matrices
the_file = 'edgelists/sociopatterns_hypertext.dat'
At = AdjMatrixSequence(the_file, directed=False, write_label_file=False)
#At.info_scipy_version()

# compute accessibility
c = At.unfold_accessibility()

# derivative of accessibility profile
h = Tools.cdf2histogram(c)

# write the results to files
Tools.dict2file(c, "path_density.txt")
Tools.dict2file(h, "path_durations.txt")


### ALTERNATIVELY: read a temporal edge list and randomize it
### Details about randomization techniques are in Supplementary Material of
### [1] Lentz et al, Phys. Rev. Lett. 110, 118701 (2013)
###

# import a temporal network edgelist for randomization
# E=TemporalEdgeList(the_file,directed=False)
# E.RE()
# E.write("Randomized_edges_RE.txt")
