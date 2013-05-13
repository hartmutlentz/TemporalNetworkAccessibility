#!/usr/bin/env python
#
#
import sys
sys.path.append('./src')
from AdjacencyMatrixSequence import AdjMatrixSequence
from TemporalNetworkEdgeList import TemporalEdgeList
import Tools


if __name__=="__main__":
    
    # import an edgelist as sequence of adjacency matrices
    the_file='edgelists/sociopatterns_hypertext.dat'
    At=AdjMatrixSequence(the_file,directed=True,write_label_file=False)
    
    # compute accessibility
    c=At.unfold_accessibility()
    h=Tools.cdf2histogram(c)
    Tools.dict2file(c,"path_density.txt")
    Tools.dict2file(h,"path_durations.txt")

    
    # import a temporal network edgelist for randomization
    #E=TemporalEdgeList(the_file,directed=True)
    #E.RE()
    #E.write("Test_RE.txt")



