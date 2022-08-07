---
layout: single
title: "CONNECTOR - How to install"
permalink: /HowInstall/
toc: true
toc_label: "Requirements"
toc_icon: "cog"
---
<style type="text/css">
  body{
  font-size: 1.125em;
}
</style>

## CONNECTOR

it is built under R 4.1.0, for previous versions we suggest installing the right version of the package *locfit*:
```
install_version("locfit", version = "1.5-9.2")
```

Note that R should be built with tcltk support. To check this, you should run *capabilities("tcltk")* from R. If it returns FALSE then we suggest checking how to cope with this from [here](https://stackoverflow.com/questions/25212800/error-onload-failed-in-loadnamespace-for-tcltk).

Before installing **CONNECTOR**, the following R packages have to be installed:

```
install.packages(c('cowplot', 'fda', 'flexclust','ggplot2', 'MASS',
                   'Matrix', 'plyr','ggplotify', 'RColorBrewer',
                   'readxl', 'reshape2', 'splines', 'statmod', 
                   'sfsmisc','shinyWidgets', 'viridis', 'dashboardthemes',
                   'shinybusy','shinydashboard','shinyjs','tidyr',
                   'shinyFiles','devtools'))"
```

Finally, to install **CONNECTOR** you can use the **devtools** R package:

```
library(devtools)
install_github("https://github.com/qBioTurin/connector", ref="master")
```

Here the CONNECTOR step by step guide: [CONNECTOR pdf](https://github.com/qBioTurin/connector/raw/master/vignettes/CONNECTORguide.pdf).

## Docker
If you want to use Docker (which is not mandatory), you need to have it installed on your machine. For more info see this document:
https://docs.docker.com/engine/installation/.

### How to install docker
Ensure your user has the rights to run docker (without the use of ```sudo```). To create the docker group and add your user:

* Create the docker group.

```
  $ sudo groupadd docker
```
* Add your user to the docker group.

```
  $ sudo usermod -aG docker $USER
```
* Log out and log back in so that your group membership is re-evaluated.

### How to run CONNECTOR from docker
* **Without installing CONNECTOR:** 

 1. Download the image from [dockerhub](https://hub.docker.com/layers/qbioturin/connector/latest/images/sha256:47a2db335ce28139952542ebe86a067373e1e7a3fafb1003b73786b52a26b96f)
```
  $ docker pull qbioturin/connector:latest
```
 2. Run the RStudio Docker image with CONNECTOR installed and work with RStudio in the browser
 ```
  $ docker run --cidfile=dockerID -e DISABLE_AUTH=true -p 8787:8787 -d qbioturin/connector:latest
  $ open http://localhost:8787/
 ```
 
* **From CONNECTOR:** 
 it is possible to download and run the RStudio Docker image by using the CONNECTOR functions called *downloadContainers* and *docker.run* as follows:
 ```
 > downloadContainers()
 > docker.run()
 ```
