FROM rocker/tidyverse:latest

RUN apt-get update -y
RUN apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev
RUN apt-get install -y curl libcurl4-openssl-dev

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

COPY . /agoradataprocessing
WORKDIR /agoradataprocessing

RUN Rscript -e 'devtools::install_deps(pkg = ".", dependencies = TRUE, threads = getOption("Ncpus",1))'
RUN R CMD INSTALL .

COPY exec/process.R /usr/local/bin/
