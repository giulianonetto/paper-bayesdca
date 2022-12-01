FROM giulianocruz/rstudio:0.0.4

RUN R -e "devtools::install_version('targets', version = '0.14.0', dependencies = T)"
RUN R -e "devtools::install_version('OptimalCutpoints', version = '1.1.5', dependencies = T)"
RUN R -e "devtools::install_github('giulianonetto/bayesdca@d6fadb8efc186ddd99e999d25479b3489948d537')"