FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@5295e3c320f747cf50839b2aa3ddae936e69f45e')"
