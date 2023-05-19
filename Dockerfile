FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@ef1778b39c6f86735cbfabe779843c087929be25')"
