FROM rocker/rstudio:4.1.0

RUN apt-get clean all && \
	apt-get update && \
	apt-get upgrade -y && \
	apt-get install -y \
		libhdf5-dev \
		libcurl4-gnutls-dev \
		libssl-dev \
		libxml2-dev \
		libpng-dev \
		libxt-dev \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev \
		libglpk40 \
		libgit2-28 \
	&& apt-get clean all && \
	apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install all required R packages
RUN Rscript -e "install.packages(c('cowplot', 'fda', 'flexclust','ggplot2', 'MASS', 'Matrix', 'plyr','ggplotify', 'RColorBrewer', 'readxl','reshape2', 'splines', 'statmod', 'sfsmisc','shinyWidgets', 'viridis', 'dashboardthemes','shinybusy','shinydashboard','shinyjs','tidyr'))"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('shinyFiles')"
# Create scratch directory
#RUN mkdir /home/scratch; chmod 777 /home/scratch

# Create data directory
#RUN mkdir /home/data; chmod 777 /home/data

# date > marker
# COPY marker /dev/null

RUN Rscript -e " devtools::install_github('qBioTurin/connector', ref='master',dependencies=TRUE, force = TRUE)"

