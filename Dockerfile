FROM giulianocruz/rstudio:0.0.9
RUN R -e "devtools::install_github('giulianonetto/bayesdca@5295e3c320f747cf50839b2aa3ddae936e69f45e')"
RUN R -e "devtools::install_version('predtools', version = '0.0.3', dependencies = T)"
RUN echo "force_color_prompt=yes" >> ~/.bashrc
ENV R_REMOTES_UPGRADE=never
ENV VIRTUAL_ENV_DISABLE_PROMPT=1
RUN echo "PS1='\[\e[1;38;2;231;41;138m\]${VIRTUAL_ENV:+[$(basename -- $VIRTUAL_ENV)] }\[\e[1;38;2;117;112;179m\][[\u]]\[\033[00m\]:\[\e[1;38;2;27;158;119m\]\w/\n\[\e[1;38;2;217;95;2m\]\\$\\$\[\033[00m\] '" >> ~/.bashrc
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    openssh-client
RUN echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections
RUN apt-get update \
    && apt-get install -y --no-install-recommends ttf-mscorefonts-installer \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
RUN fc-cache -fv