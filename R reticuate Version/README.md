TemporalNetworkAccessibility
================
-   [Overview](#overview)
-   [Package installation and Code download](#package-installation-and-code-download)
-   [Run code](#Run-code)



# Overview

This code shows how to run the Python code of the  [TemporalNetworkAccessibility](https://github.com/hartmutlentz/TemporalNetworkAccessibility) in R using the package *reticulate* .

It provides code for computing the results of the Papers

**Unfolding Accessibility Provides a Macroscopic Approach to Temporal Networks**,  
Lentz et al., Phys. Rev. Lett. 110, 118701 (2013)

and

**Disease Spread through Animal Movements: A Static and Temporal Network Analysis of Pig Trade in Germany**,
Lentz et al., PLOS ONE 11, e0155196–32 (2016)

Please cite the second reference (PLOS ONE) if you use this R code.


A full description and the *Python* code can be found at [TemporalNetworkAccessibility](https://github.com/hartmutlentz/TemporalNetworkAccessibility)

Additionally, we provide a tutorial similar to the tutorial for the Python code.


# Package installation and Code download

## Python

You need an installed version of [Python 3](https://www.python.org/downloads/).
In Windows, I had best experiences installing it in a separate forlder (e.g. E://Python3/).
Please install the libraries *pandas*,*numpy*,*networkx*, *scipy*, and *matplotlib*

## Code download

Download the Python code from [TemporalNetworkAccessibility](https://github.com/hartmutlentz/TemporalNetworkAccessibility).
Go to the main folder of the extracted download.

## R Packages
To run the functions you need to install the following libraries 

``` r
install.packages("reticulate", dependencies = TRUE)
``` 
Then load the packages ```reticulate``` and check, if python is found.
Eventually define the folder, where Python is installed.

``` r
library(reticulate)
## check if python is available
Sys.which("python")
use_python("E://Python3")

``` 

## Python Libraries

``` r
## load python packages
os <- import("os")
pd <- import("pandas")
scipy <- import("scipy")
sp <- import("scipy.sparse")
np <- import("numpy")
plt <- import("matplotlib")
plt <- import("matplotlib.pyplot")
nx <- import("networkx")

setwd("./src")
source_python("AdjacencyMatrixSequence.py")
source_python("Tools.py")
source_python("TemporalNetworkEdgeList.py")
setwd("../")
``` 

# Run code
## Load data

``` r
At = AdjMatrixSequence("./edgelists/sociopatterns_hypertext.dat","directed")
```

## Compute Accessibility and its derivative

```{r, message=FALSE, warning=FALSE}

c = AdjMatrixSequence$unfold_accessibility(At)
#c1 = AdjMatrixSequence$unfold_accessibility_memory_efficient(At)

h = np$gradient(c)

```

## Plot the results
```{r,  message=FALSE, warning=FALSE}

plot(c, type = "l", col="blue",
     xlab = "time", 
     ylab = "cumulative #paths")


plot(h, type = "l",
     xlab = "time", 
     ylab = "# shortest paths")

```

## Compute causal fidelity

```{r causal fidelity,  message=FALSE, warning=FALSE}
# Causal fidelity
causal_paths = c[length(c)]
static_paths = At$static_path_density()
print(paste0("---> Causal fidelity is ", (causal_paths)/(static_paths)))

```


## Compute Causal Fidelity for an increasing aggregation window 
Path density of the static network for different aggregation windows
Warning: this might be slow. In doubt set verbose=True to see the progress.

```{r causal fidelity1,  message=FALSE, warning=FALSE}
c2 = AdjMatrixSequence$step_by_step_static_path_density(At, verbose="False")

# causal fidelity is number of causal paths divided by number of static paths
c_ff <- vector()
for (i in 1:length(c)){
  c_ff[i] = c[i]/c2[i]}

# plot the result

plot(h, type = "l",
     xlab = "aggregation window [time]", 
     ylab = "causal fidelity")

```
Considering the network as static starts to make sense after about 2000 time steps. The high causal fidelity at the beginning is a trvial effect: the network is very sparse and most paths are short (length 1), so causality does not play a role yet.


## Tracing forward and backward
This is just single node accessibility unfolding. We use pandas here for plotting, and matplotlib works, too.

```{r Tracing,  message=FALSE, warning=FALSE}

##### tracing forward ####
node_name <- 101
c1 = AdjMatrixSequence$unfold_accessibility_single_node(At, node_name)
c1df = pd$DataFrame(c1)


plot(c1, type = "l", col = "blue",
     xlab = "time",
     ylab = "affected nodes")

# tracing backward
Bt = At
Bt$time_reversed()
c2 = AdjMatrixSequence$unfold_accessibility_single_node(Bt, node_name)
c2df = pd$DataFrame(c2)
plot(c2, type = "l", col = "darkblue",
     xlab = "time",
     ylab = "affected nodes")

```


if you want to have the nodes reachable from the starting node, use `At.trace_forward(node_name)`. Note that the resulting dictionary can be huge.

```{r Tracing1,  message=FALSE, warning=FALSE}

node_name = 10
stop_node =1000
reachable_nodes = AdjMatrixSequence$trace_forward(At, start_node = as.integer(node_name), stop = as.integer(stop_node))
print(paste0("The reachable nodes of node ",node_name," after ",stop_node,"time steps are:",reachable_nodes[[stop_node]]))
```

## Unfold Accessibility

```{r}

the_file = "./edgelists/sexual_contacts.dat"
At = AdjMatrixSequence(the_file, directed="True", write_label_file="False")


c = AdjMatrixSequence$si_model(At, p = 0.1)
c01 = AdjMatrixSequence$si_model(At, p = 0.01)
c99 = AdjMatrixSequence$si_model(At, p = 0.99)


plot(c99, type = "l", col = "black",
     xlab = "time",
     ylab = "affected nodes")
lines(c, col = "blue")
lines(c01, col = "green")

legend("topleft", c("p=0.99", "p=0.1","p=0.01"),
       lty = c(1,1),
       col = c("black", "blue", "green"))
```


