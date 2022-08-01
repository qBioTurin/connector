---
layout: single
title: "CONNECTOR - How to install"
permalink: /HowInstall/
toc: true
toc_label: "Requirements"
toc_icon: "cog"
---


## CONNECTOR
To install **CONNECTOR** you can use the **devtools** R package:

```
install.packages("devtools")
library(devtools)
install_github("https://github.com/qBioTurin/connector", ref="master")
```

Here the CONNECTOR step by step guide: [CONNECTOR pdf](https://github.com/qBioTurin/connector/raw/master/vignettes/CONNECTORguide.pdf).

## Docker
If you want to use Docker (which is not mandatory), you need to have it installed on your machine. For more info see this document:
https://docs.docker.com/engine/installation/.

Ensure your user has the rights to run docker (witout the use of ```sudo```). To create the docker group and add your user:

* Create the docker group.

```
  $ sudo groupadd docker
```
* Add your user to the docker group.

```
  $ sudo usermod -aG docker $USER
```
* Log out and log back in so that your group membership is re-evaluated.

