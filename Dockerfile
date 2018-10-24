FROM rocker/tidyverse:3.4.0

RUN R -e "devtools::install_github('Sage-Bionetworks/agoradataprocessing')"
