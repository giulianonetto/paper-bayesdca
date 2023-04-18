FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@4dd2c7840f61f1d339ee765b1f735dd491a5c852')"
