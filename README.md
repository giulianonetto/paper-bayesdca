# paper-bayesdca
Analyses for the bayesDCA paper

To run the pipeline step named `"sample_size_simulation"`, run:

```
docker build -t bayesdcapaper . && \
    docker run -it --rm -v ${PWD}:/home/rstudio bayesdcapaper \
    R -e "targets::tar_make(names = 'sample_size_simulation')"
```

or just remove `names = 'sample_size_simulation'` to run all steps (might take a while).