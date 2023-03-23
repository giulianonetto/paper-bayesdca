FROM giulianocruz/rstudio:0.0.8
RUN R -e "devtools::install_github('giulianonetto/bayesdca@48e31e8414d8fca64f9447c1c6328b1c931b623d')"