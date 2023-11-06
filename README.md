# connector
It is R package that is able to fit and cluster growth data using four different fitting models (i.e. Functional Clustering Models, Malthus model, Gompertz model, and logistic model). 
It exploits an unsupervised clustering algorithm to cluster the fitted data, and the separation and tightness measures are provide to evaluete the quality of the derived  clusters.

## Required installed packages
The following R packages must be installed to use connector:
cowplot, fda, flexclust, ggplot2, MASS, Matrix, plyr, ggplotify, RColorBrewer, readxl, reshape2, splines, statmod, sfsmisc, shinyWidgets, viridis and dashboardthemes.

```
install.packages(c("cowplot", "fda", "flexclust", "ggplot2", "MASS", "Matrix", "plyr", "ggplotify",
	"RColorBrewer", "readxl", "reshape2", "splines", "statmod", "sfsmisc", "shinyWidgets", "viridis", "dashboardthemes"))
```

## How to install connector
To install it you can use  **devtools**:

```
install.packages("devtools")
library(devtools)
install_github("qBioTurin/connector", ref="Classification",dependencies=TRUE)
```

## An example of connector analysis:
The **MainCommandsList.R** script contains the list of commands that can be used to reproduce the analysis described in the paper "CONNECTOR: fitting and clustering analysis of biological growth data".

To execute this script you can use **demo()**
```
library(connector)
demo("MainCommandsList", package = "connector")
```
Observe that some functions in the **MainCommandsList.R** script require several minutes to be terminated. 

## Diclaimer:
CONNECTOR developers have no liability for any use of docker4seq functions, including without limitation, any loss of data, incorrect results, or any costs, liabilities, or damages that result from the use of CONNECTOR. 

## How to cite

```
@article{pernice2023connector,
  title={CONNECTOR, fitting and clustering of longitudinal data to reveal a new risk stratification system},
  author={Pernice, Simone and Sirovich, Roberta and Grassi, Elena and Viviani, Marco and Ferri, Martina and Sassi, Francesco and Alessandr{\`\i}, Luca and Tortarolo, Dora and Calogero, Raffaele A and Trusolino, Livio and others},
  journal={Bioinformatics},
  volume={39},
  number={5},
  pages={btad201},
  year={2023},
  publisher={Oxford University Press}
}
```
