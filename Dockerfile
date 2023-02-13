FROM giulianocruz/rstudio:0.0.7
RUN R -e "devtools::install_github('giulianonetto/bayesdca@fcb7b7e152aff220ae3c89ba329ec0d3565008ec')"