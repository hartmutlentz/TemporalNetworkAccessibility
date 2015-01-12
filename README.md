# TemporalNetworkAccessibility

Provides classes used for computing the results of the Paper  
**Unfolding Accessibility Provides a Macroscopic Approach to Temporal Networks**, Lentz et al., Phys. Rev. Lett., 2013.

## Using the Software
You can compute the Accessibility Matrix of a temporal netwotk like this:
### Step 1
```python
A = AdjMatrixSequence(<your input file>)
```
This generates a list of adjacency matrices, i.e. a temporal network in matrices representation.

### Step 2
```python
c = A.unfold_accessibility()
```
This computes eqn. (4) in the Paper, that is $$ \mathcal{P}_n = \bigwedge _i (\mathbf{1} \vee \mathbf{A}_i) $$ step by step and returns a dictionary containing the path density (as black line in Fig. 2).


An example of usage is shown in 'Unfold_Accessibility.py'.




Required Software/packages:
- Python 2.7 (should run on 2.6, too)
- scipy
- numpy
- Networkx package (optional, required for the configuration model).

The class AdjMatrixSequence is implemented using scipy.sparse matrices.
These matrices are restricted to 2^31 nonzero entries, regardless of the used memory!




Datasets:
Exemplary temporal network datasets are provided in the 'edgelists' folder.
The file 'edgelists/sexual_contacts.dat' is from [1] and the file 'edgelists/sociopatterns_hypertext.dat' can be downloaded from sociopatterns.org [2].
Both files were also used in [3].

[1]	L. E. C. Rocha, F. Liljeros, and P. Holme, Proc. Natl. Acad. Sci. U.S.a. 107, 5706 (2010).
[2]	L. Isella, J. Stehle, A. Barrat, C. Cattuto, J.-F. Pinton, and W. Van den Broeck, J. Theor. Biol. 271, 166 (2011).
[3]	H. H. K. Lentz, T. Selhorst, and I. M. Sokolov, Phys. Rev. Lett. 110, 118701 (2013).

