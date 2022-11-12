FROM giulianocruz/rstudio:0.0.4

RUN R -e "devtools::install_version('targets', version = '0.14.0', dependencies = T)"