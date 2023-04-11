FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@50097b2d66636f0240ce4ebfea9e3b44e0926f0f')"