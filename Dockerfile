FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@715eb809c5cb041a52b607a6a588b25302e31db0')"