FROM python:3.8

RUN apt-get update && apt-get install -y --no-install-recommends r-base build-essential
RUN apt-get upgrade libstdc++6 -y
RUN add-apt-repository pa:ubuntu-toolchain-r/test
RUN apt-get update 
RUN apt-get upgrade
RUn apt-get dist-upgrade
RUN R -e "install.packages(c('argparse', 'R.utils', 'data.table', 'ggplot2', 'stringr', 'reshape2', 'dplyr', 'plyr', 'ggpubr'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "update.packages(ask = FALSE, checkBuilt = TRUE)"