FROM giulianocruz/rstudio:0.0.8
RUN R -e "devtools::install_version('simsurv', version = '1.0.0', dependencies = FALSE, upgrade_dependencies = FALSE)"