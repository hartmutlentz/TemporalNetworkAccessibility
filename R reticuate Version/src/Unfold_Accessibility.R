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

## load python scripts
setwd("./src")
source_python("AdjacencyMatrixSequence.py")
source_python("Tools.py")
source_python("TemporalNetworkEdgeList.py")
setwd("../")

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

DF <- data.frame(p99 = AdjMatrixSequence$si_model(At, p = 0.99),
                 p10 = AdjMatrixSequence$si_model(At, p = 0.1),
                 p01 = AdjMatrixSequence$si_model(At, p = 0.01))
