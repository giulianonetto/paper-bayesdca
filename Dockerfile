FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@9056cf356dad2fd06d0f9e244b9c7f114255978a')"
