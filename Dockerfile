FROM giulianocruz/rstudio:0.0.4

RUN R -e "devtools::install_version('targets', version = '0.14.0', dependencies = T)"
RUN R -e "devtools::install_version('OptimalCutpoints', version = '1.1.5', dependencies = T)"
RUN R -e "devtools::install_github('giulianonetto/bayesdca@cf53183d22aba7b54bcd8c664b68ca83973700a2')"