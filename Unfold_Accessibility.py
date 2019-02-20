#!/usr/bin/env python
#
#
import sys
import numpy as np
sys.path.append('./src')
from AdjacencyMatrixSequence import AdjMatrixSequence
from TemporalNetworkEdgeList import TemporalEdgeList
import Tools

# import an edgelist as sequence of adjacency matrices
the_file = 'edgelists/Test.dat'
the_file = "edgelists/sexual_contacts.dat"
At = AdjMatrixSequence(the_file, directed=False, write_label_file=False)

# compute accessibility
c = At.unfold_accessibility(return_accessibility_matrix=False)

# derivative of accessibility profile
h = np.gradient(c)

# write the results to files
Tools.dict2file(c, "shortest_path_durations_cumulative.txt")
Tools.dict2file(h, "shortest_path_durations_histogram.txt")

# Causal fidelity
causal_paths = c[-1]
static_paths = At.static_path_density()
print("---> Causal fidelity is ", float(causal_paths)/float(static_paths))

### ALTERNATIVELY: read a temporal edge list and randomize it
### Details about randomization techniques are in Supplementary Material of
### [1] Lentz et al, Phys. Rev. Lett. 110, 118701 (2013)
###

# import a temporal network edgelist for randomization
#E = TemporalEdgeList(the_file, directed=True)
#E.RE(5)
#E.write("Randomized_RE.txt")
