FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@810977adbb49a9be99e7adfdadafb341cb119589')"
