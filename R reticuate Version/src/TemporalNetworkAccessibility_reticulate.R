## ----load packages, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------------------
#Sys.setenv(RETICULATE_PYTHON = "E://Python3/python")
library(reticulate)

## check if python is available
Sys.which("python")
use_python("E://Python3")


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
#source_python("Unfold_Accessibility.py")



## ---- message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------------------------
At = AdjMatrixSequence("./edgelists/sociopatterns_hypertext.dat","directed")



## ---- message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------------------------

c = AdjMatrixSequence$unfold_accessibility(At)
#c1 = AdjMatrixSequence$unfold_accessibility_memory_efficient(At)

h = np$gradient(c)



## ----  message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------------------

plot(c, type = "l", col="blue",
     xlab = "time", 
     ylab = "cumulative #paths")


plot(h, type = "l",
     xlab = "time", 
     ylab = "# shortest paths")



## ----causal fidelity,  message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------------
# Causal fidelity
causal_paths = c[length(c)]
static_paths = At$static_path_density()
print(paste0("---> Causal fidelity is ", (causal_paths)/(static_paths)))



## ----causal fidelity1,  message=FALSE, warning=FALSE------------------------------------------------------------------------------------------------------------
c2 = AdjMatrixSequence$step_by_step_static_path_density(At, verbose="False")

# causal fidelity is number of causal paths divided by number of static paths
c_ff <- vector()
for (i in 1:length(c)){
  c_ff[i] = c[i]/c2[i]}

# plot the result

plot(h, type = "l",
     xlab = "aggregation window [time]", 
     ylab = "causal fidelity")



## ----Tracing,  message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------

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



## ----Tracing1,  message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------

node_name = 10
stop_node =1000
reachable_nodes = AdjMatrixSequence$trace_forward(At, start_node = as.integer(node_name), stop = as.integer(stop_node))
print(paste0("The reachable nodes of node ",node_name," after ",stop_node,"time steps are:",reachable_nodes[[stop_node]]))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------

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

