# TemporalNetworkAccessibility

Provides classes used for computing the results of the Paper  
**Unfolding Accessibility Provides a Macroscopic Approach to Temporal Networks**,  
Lentz et al., Phys. Rev. Lett. 110, 118701 (2013).

Please cite this reference, if you use the software.

This software has been written for **Python 2** and discontinued. If you are using Python 3 (recommended), please choose the branch ```master``` on this website.

## Using the Software
To compute the results of the paper, you only need the class *AdjMatrixSequence*.
You can compute the Accessibility Matrix of a temporal network following these steps:

### Step 1
```python
A = AdjMatrixSequence("<your input file>")
```
This generates a list of adjacency matrices, i.e. a temporal network in matrices representation. The input file ```<your input file>``` is a simple text file with 3 columns and must be in the following form:
```
source_node	target_node	time
```
where the columns are separated by tabs. Thus, the file is a standard edge list with an additional time column. Examples are in the *edgelists* folder.

Additionally, more options could be used in this step. You can use option *directed=True*, if your network is directed (default is undirected). The option *write_label_file=True* (default is *False*) can be used, if you have nodes names that are not "matrix label friendly", as names for example. If *write_label_file=True*, an additional file is written into the working directory, which contains the old and new node labels.

### Step 2
```python
c = A.unfold_accessibility()
```
This computes eqn. (4) in the Paper, that is  
$ \mathcal{P}_n = \bigwedge _i (\mathbf{1} \vee \mathbf{A}_i)$,  
step by step and returns a dictionary containing the path density. In the paper [3] this is the black line in Fig. 2 for example.

**Memory efficience**. When using ```A.unfold_accessibility()``` like above, computations are done using matrices. These matrices can become huge for large networks. To save memory, you can use the following variant:
```python
c = A.unfold_accessibility_memory_efficient()
```
This variant can be slower, but works also for large networks. It decomposes the matrices into vectors and computes the solution for one vector at a time. (This method can be parallelized (to be done)).

### Step 3
```python
h = np.gradient(c)
```
This returns the numerical derivative of the path-density, which is the shortest-path-duration-distribution. In the paper [3] this is the red line in Fig. 2 for example.

### Step 4
```python
Tools.dict2file(c, "path_density.txt")
Tools.dict2file(h, "path_durations.txt")
```
This writes the generated data to txt-files, so you can plot it using gnuplot, Excel or you favorite plotting software. The files have 2 columns: *time and path density* or *time and path duration*, respectively.

A working example with steps 1—4 is shown in the file 'Unfold_Accessibility.py'.

**A note on normalization.** Path-density (```c```) and the distribution of shortest path durations (```h```) are not normalized at this point. Consequently, you should normalize ```c``` by the number of nodes squared and ```h``` to unity. As noted in the paper, ```c``` is not necessarily normalized to unity.

### Causal fidelity
Causal fidelity is given by the causal path density normalized by the static path density. The causal path density is simply given by the last computed value ```c[-1]```. The static path density can be computed using method ```At.static_path_density()```.

After ```c``` has been computed, the code would look like this:

```python
causal_paths = c[-1]
static_paths = At.static_path_density()
print "---> Causal fidelity is ", float(causal_paths)/float(static_paths)
```

### Additional functionality
The Class *TemporalEdgeList* provides methods to load a temporal network as a temporal edgelist. It can be used for randomization of temporal networks. There is a number of methods to randomize temporal networks. The methods implemented here have also been used in the supplementary material of [3] (also see references therein). The methods described in more detail described in my PhD thesis [4].

You can create a temporal edgelist like this:
```python
E = TemporalEdgeList("<your input file>")
```
You can add the option *directed=True*, if your network is directed.
If you want to randomize the network, for example using the *randomized edges* (RE) model, you simply use
```python
E.RE()
```
Finally, you can write the new edgelist into a textfile like this:
```python
E.write("Randomized_edges_RE.txt")
```

## Required Software/packages:
- Python 2.7 (should run on 2.6 as well)
- scipy
- numpy
- Networkx package (optional, required for the configuration model).

The class AdjMatrixSequence is implemented using scipy.sparse matrices.
In older Scipy versions (< 0.14.0), sparse matrices are restricted to 2^31 (~ 10^9) nonzero entries, regardless of the used memory. Make sure, your scipy version is up to date, if you plan to handle large networks (~ 10^5 nodes.)

## Datasets
Exemplary temporal network datasets are provided in the 'edgelists' folder.
The file 'edgelists/sexual_contacts.dat' is from [1] and the file 'edgelists/sociopatterns_hypertext.dat' can be downloaded from sociopatterns.org [2].
Both files were also used in [3].

## Literature
[1]	L. E. C. Rocha, F. Liljeros, and P. Holme, Proc. Natl. Acad. Sci. U.S.a. 107, 5706 (2010).  
[2]	L. Isella, J. Stehle, A. Barrat, C. Cattuto, J.-F. Pinton, and W. Van den Broeck, J. Theor. Biol. 271, 166 (2011).  
[3]	H. H. K. Lentz, T. Selhorst, and I. M. Sokolov, Phys. Rev. Lett. 110, 118701 (2013).  
[4] H. H. K. Lentz, PhD Thesis, [Humboldt-University of Berlin](http://edoc.hu-berlin.de/dissertationen/lentz-hartmut-2013-11-06/METADATA/abstract.php?id=40377) or [GitHub](https://github.com/hartmutlentz/Thesis)
[5]	Lentz, H. H. K. et al., PLOS ONE 11, e0155196–32 (2016).
