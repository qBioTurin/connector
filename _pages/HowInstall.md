---
layout: single
title: "CONNECTOR - How to install"
permalink: /HowInstall/
toc: true
toc_label: "Requirements"
toc_icon: "cog"
---


## CONNECTOR
To install **CONNECTOR** you can use use **devtools**:

```
install.packages("devtools")
library(devtools)
install_github("https://github.com/qBioTurin/connector", ref="master")
```


## Docker

You need to have docker installed on your machine, for more info see this document:
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

