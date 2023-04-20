FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@7fb22fe9d6dfdcccebb16dcb7122e10a6d212f18')"
