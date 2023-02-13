---
title: "Applied Survival Models"
author: "Jacqueline Buros Novik"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Applied survival models}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, eval = T, results = 'hide', echo = F}
knitr::opts_chunk$set(message = FALSE,
                      warning = FALSE,
                      fig.width = 6,
                      fig.height = 6
                      )
```


```{r load-packages, eval = T, echo = F}
library(httr)
library(readr)
library(cgdsr)
library(purrr)
library(dplyr)
library(assertthat)
library(ggplot2)
library(survival)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))
library(shinystan)
library(gridExtra)
library(ggfortify)
library(scales)
library(stringr)
```

Survival modeling is a core component of any clinical data analysis toolset.

Here we will work through an example of fitting a survival model in Stan, using as an example data from TCGA on patients with Bladder Urothelial Carcinoma. 

Fitting survival models in Stan is fairly straightforward. However, specifying a model in Stan requires a bit more thought than when using standard MLE tools such as `survfit`. 

This vignette is part 1 of a multi-part series. 

# Outline

Over the course of several vignettes, we will review various models and approaches to analysis. 

In this vignette, we cover the following:

1. Data review & inspection
2. Review parametric Weibull survival model 
3. Test parametric survival model against simulated data
4. Fit NULL parametric survival model with TCGA data
5. Check convergence & review posterior predictive checks on model

Part II of this series will consider a non-parametric piecewise exponential model.

# Data

We will use data from [The Cancer Genome Atlas (TCGA)](https://tcga-data.nci.nih.gov) for patients diagnosed with Bladder Urothelial Carcinoma. This cohort has an acronym **BLCA**. 

TCGA is a repository giving access to clinical and molecular data for a variety of cancer types, primarily focusing on genomic datasets and outcomes following standard therapeutic interventions. The BLCA cohort is described in this [this paper]().

The complete clinical & molecular data are available from TCGA data portals, but in this case we will use a curated version of these data available from <http://www.cbioportal.org/>.

To load the data into R, we *could* query these data via the web-url service: 

```{r example-load-data, eval = FALSE}
url <- 'http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id=blca_tcga_all'
req <- httr::GET(url)
clinical_data <- 
    httr::content(req,
                  type = 'text/tab-separated-values',
                  col_names = T,
                  col_types = NULL
                  )
str(clinical_data)
```

However, MSKCC has provided the [CGDS-R](http://cran.r-project.org/web/packages/cgdsr/index.html) package, which provides an easier interface to the same data.

Run the following lines of code to load the clinical data for this cohort:

```{r actual-load-data}
mycgds = cgdsr::CGDS("http://www.cbioportal.org/public-portal/")
selected_case_list = 'blca_tcga_all'
clinical_data = cgdsr::getClinicalData(mycgds, selected_case_list)
```

Describe the structure of data available:

```{r inspect-data}
str(clinical_data,  no.list = T, vec.len = 2)
```

## Data cleaning 

Let's do some minimal data manipulation on this DataFrame to make it easier to work with. I personally happen to like lower-case variable names and to work with NA values instead of empty strings.

We will save the converted data in a data.frame called `clin_data`.

```{r initprep-data}
## names to lower case
names(clinical_data) <- tolower(names(clinical_data))

## convert empty strings -> NA values
convert_blank_to_na <- function(x) {
    if (!purrr::is_character(x)) {
        warning('input vector is not character - returning original input')
        return(x)
    } else {
        ifelse(x == '', NA, x)
    }
}
clin_data <- clinical_data %>%
    dplyr::mutate_each(funs = funs(convert_blank_to_na), everything())

## inspect resulting data frame
str(clin_data, vec.len = 2, list.len = 10)
```

## Data exploration

In survival analysis, the outcome or dependent variable is the *time to event* where some event times are not observed (IE they are censored). 

Here we consider the more common scenario of *right-censoring*. This is the case where the terminating event is not observed. Observations are instead censored at time `t` .

Our first analysis will treat **overall survival** as the event of interest, as opposed to progression-free survival. In this cohort, the overall survival is described by two variables: `os_status` & `os_months`. 

We will start by inspecting these data.

### Negative/missing event times

We have one observation with missing data for `os_status`:

```{r inspect-os-status}
clinical_data %>%
    dplyr::filter(is.na(os_status) | os_status == '') %>%
    dplyr::select(os_status, os_months) %>%
    str()
```

We have a 3 observations with unknown and/or negative survival times (`os_months < 0`):

```{r inspect-os-months}
clinical_data %>%
    dplyr::filter(!is.na(os_status) & os_status != '') %>%
    dplyr::filter(os_months < 0 | is.na(os_months)) %>%
    dplyr::select(os_status, os_months) %>%
    head()
```

For now, we remove these observations from our analysis dataset (`clin_data`).

```{r create-clin-data}
clin_data <- 
    clin_data %>%
    dplyr::filter(!is.na(os_status) & os_status != '') %>%
    dplyr::filter(os_months >= 0 & !is.na(os_months))

## confirm 4 fewer observations than original
assert_that(nrow(clin_data) == nrow(clinical_data) - 4)
```

### Distribution of event times

Among the remaining observations, the median time to event is `r median(clin_data$os_months)` months.

Event times are distributed as follows among observed (DECEASED) & censored (LIVING) observations:

```{r plot-os-months}
ggplot(clin_data,
       aes(x = os_months,
           group = os_status,
           colour = os_status,
           fill = os_status
           )) + 
    geom_density(alpha = 0.5)
```

The KM curve for these observations additionally looks like this: 

( constructed using `survfit` )

```{r km-survival-curve}
mle.surv <- 
    survfit(
        Surv(os_months,os_deceased) ~ 1,
        data = clin_data %>%
            dplyr::mutate(os_deceased = os_status == 'DECEASED')
    )
autoplot(mle.surv, conf.int = F) +
    ggtitle('KM survival curve for BLCA cohort')
```


# First analysis: parametric survival model

For our first analysis we will work with a parametric Weibull survival model. 

We will start with model code adapted from  [wei_bg.stan](https://github.com/to-mi/stan-survival-shrinkage/blob/master/wei_bg.stan) within the [github repo]('http://github.com/to-mi/stan-survival-shrinkage') accompanying [Peltola et al, 2014](http://ceur-ws.org/Vol-1218/bmaw2014_paper_8.pdf)'s nice paper describing a bayesian approach to biomarker evaluation.

This model assumes that the time to event `x` follows a Weibull distribution. 

Stan parameterizes this probability density function as :

$$f(x|\alpha,\sigma) = 
\frac{\alpha}{\sigma}\left(\frac{x}{\sigma}\right)^{\alpha-1}e^{-(x/\sigma)^{\alpha}}$$

In the context of this analysis, we will define two parameters:

* `alpha` (shape) defined as above
* `mu` (scale) where $\sigma = e^\frac{-\mu}{\alpha}$.

If we had covariates and wanted to estimate a proportional hazards model, we would replace `mu` with a linear combination of covariates. However, in this case we are interested in recovering features of our NULL model and so we treat `mu` as a constant intercept.

## Stan code for the model

The stan code for this model is included in this [biostan package](http://github.com/jburos/biostan) as `weibull_survival_null_model.stan`. 

It can be accessed by calling `system.file()`, as:

```{r locate-stan-file}
if (!require(biostan))
    devtools::install_github('jburos/biostan')
library(biostan)
stan_file <- system.file('stan', 'weibull_survival_null_model.stan', package =  'biostan')
```

Here are the contents of this file:

```{r print-stan-code}
biostan::print_stan_file(stan_file)
```

### The model in detail

Before using this model for analysis, we want to first review the model code in detail & test it against some simulated data. 

This will ensure that (a) we understand the model well, and (b) the model can recover estimates from simulated data. 

*( As you will see, several parts of the simulate-data process can also be re-used for posterior predictive checking. So we will save components of the process to be reused in later steps. )*

If you're at an R console, you can open the Stan file in an editor as follows:

```{r edit-stan-file, eval = F}
if (interactive())
    file.edit(stan_file)
```

#### Review data block 

Let's review the data block first. 

This will tell us the structure and format of data input to the model.

```{r view-data-block}
print_stan_file(stan_file, section = 'data')
```

The censored & observed data points are provided as separate input vectors. 

*observed data points*

* `Nobs`: number of observed data points 
* `yobs`: times to observed events

*censored data points*

* `Ncen`: number of censored data points 
* `ycen`: times to censored events

Recall that this is a NULL model (with no covariate values), so the number & values of observed covariates are not needed.

#### Review model block

The stan code contains an implicit constant term, in the linear predictor `mu`.

```{r print-model-block}
print_stan_file(stan_file, section = 'model')
```

Observe how the ccdf (complementary cumulative distribution function) is used 
to compute the log probability of the censored observations. 

*What does the ccdf represent in this scenario?*

*How does the model address the censoring process?*

#### Review parameters block

Our stan code also contains a reparameterization of the `alpha` term, in the `transformed parameters` block. 

Observe:

```{r print-parameters-block}
print_stan_file(stan_file, section = 'transformed parameters')
```

(recall that `tau_al` is a constant scaling term set to 10, and `alpha_raw` is a parameter with a normal(0, 1) prior distribution).

This reparameterization achieves two things : 

1. The use of `tau_al * alpha_raw` is an example of a **non-centered parameterization**. 
    - It would have been *mathematically* equivalent to define a (non-transformed) parameter `alpha` with a prior `normal(0, 10)`. 
    - However, this parameterization yields a parameter (`alpha_raw`) which is on a similar scale as other parameters in our model. The `exp()` transformation makes the difference between these two scales even more dramatic.
    - In general, having all parameters on a similar scale makes the sampling more efficient.


2. The `exp()` transformation of this parameter additionally allows us to put a prior on `log alpha`. 
    - we want to put a prior on `log alpha` since alpha enters into our likelihood in the exponent.

This seems like a lot of gymnastics to be doing. 

However, it has practical implications for our modeling efficiency.

Observe that, for a single value of `alpha_raw` (e.g. 0.2), the transformation yields:

```{r single-value-alpha}
alpha_raw <- 0.2
tau_al <- 10
log_alpha <- alpha_raw * tau_al
alpha <- exp(log_alpha)
print(alpha)
```

which may seem silly.

**However**

Consider the resulting distribution of alpha over a range of values for `alpha_raw` sampled from our `normal(0, 1)` prior:

```{r dist-alpha}
alpha_raw <- rnorm(1000, 0, 1)
tau_al <- 10
log_alpha <- alpha_raw * tau_al
alpha <- exp(log_alpha)
ggplot(data.frame(alpha = alpha, alpha_raw = alpha_raw), 
       aes(x = alpha)) + 
    geom_density() + 
    scale_x_log10(labels = scientific)
```

Notice how `alpha` ranges from 1e-10 to 1e+10 on a log scale. We have to truncate this dramatically to even consider plotting it on its original scale. 

Sampling this parameter space may require different step sizes & different tuning parameter values throughout this distribution.

The `alpha_raw` scale, by comparison, is a lot friendlier. 

```{r dist-alpha-vs-raw}
ggplot(data.frame(alpha = alpha, alpha_raw = alpha_raw), 
       aes(x = alpha, y = alpha_raw)) + 
    geom_density2d() + 
    scale_x_log10(labels = scientific)
```

This distribution is centered at 0 and has more consistent behavior throughout its range of values.

What's important to note here is that while the non-centered parameterization is mathematically equivalent to the standard parameterization, it is (in some ways) a *different model*. Consider that the reparameterization will impact the posterior estimates of parameter values.

Packages like [rstanarm](http://github.com/stan-dev/rstanarm) which provide easy wrappers to a variety of standard models implemented in Stan use a non-centered parameterization by default.

More information on non-centered parameterization:

1. [discussion on stan-dev list](https://groups.google.com/forum/#!topic/stan-dev/9ZvhKpXlwuI)
2. [Gelman, 2004. Parameterization and Bayesian Modeling](http://www.stat.columbia.edu/~gelman/research/published/parameterization.pdf)

## Testing the model on simulated data

Now that we have reviewed the model code, we are ready to simulate data according to this model.

We can simulate data using R or in Stan. We will start by simulating data in R.

### Simulate data in R 

Like our stan model code, we originally based this function on that used by the [example.R](https://github.com/to-mi/stan-survival-shrinkage/blob/master/example.R) file from the [stan-survival-shrinkage github repo](https://github.com/to-mi/stan-survival-shrinkage). 

However, after further inspection (see the related [weibull-survival-model](weibull-survival-model.html) vignette) we modified the simulate-data function slightly. Here we will work with the modified function.

Our `weibull_sim_data` function takes two parameters (`alpha` and `mu`) as inputs 
and a desired sample size (`n`). It returns a data frame of simulated event times.

```{r sim-data-function}
weibull_sim_data <- function(alpha, mu, n) {
    
    data <- data.frame(surv_months = rweibull(n = n, alpha, exp(-(mu)/alpha)),
                       censor_months = rexp(n = n, rate = 1/100),
                       stringsAsFactors = F
                       ) %>%
        dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                          'DECEASED', 'LIVING'
                                          ),
                       os_months = ifelse(surv_months < censor_months,
                                          surv_months, censor_months
                                          )
                       )

    return(data)
}
```

A few comments about this function:

1. Notice how the censoring process is `rexp()`. We chose this somewhat arbitrarily.
    - In general, our Stan model is ignorant of
the censoring process except to assume that censoring is noninformative. 
2. We have also deliberately written this function to mimic the structure of our clinical data.

This will make it easier to reuse this & other functions later.

#### Simulate data for arbitrary input values

We can use this function to simulate a dataset for hypothetical parameter values of `alpha` & `mu`.

```{r simulated-data}
test_alpha <- 0.8
test_mu <- -3

## sample size from TCGA blca data
test_n <- nrow(clin_data)

## test these inputs for arbitrary values of alpha & mu
simulated_data <- 
    weibull_sim_data(alpha = test_alpha,
                 mu = test_mu,
                 n = test_n
                 ) 
head(simulated_data)
```

The simulated survival curve looks like:

```{r sim-km-curve}
## plot KM curve from simulated data
simulated_data <- 
    simulated_data %>%
    dplyr::mutate(os_deceased = os_status == 'DECEASED')

autoplot(survival::survfit(Surv(os_months, os_deceased) ~ 1,
                      data = simulated_data
                      ), conf.int = F) + 
    ggtitle('Simulated KM curve')
```

### fit to simulated data in stan

Now that we have simulated data, we are ready to fit the model in Stan. 

If we have written both our stan code & simulated data process correctly, our posterior intervals for `alpha` and `mu` should contain the values used to simulate our dataset (`r test_alpha` and `r test_mu`).

#### preparing data for stan

Stan takes data input as a list. The contents of the list should match 
those of the `data` block in the stan code.

E.g. looking at the data block - 

```{r review-data}
print_stan_file(stan_file, section = 'data')
```

our input list to Stan should contain dimensions & values 
for observed & censored data, separately.

```{r stan-data}
observed_data <- simulated_data %>%
    dplyr::filter(os_status == 'DECEASED')

censored_data <- simulated_data %>%
    dplyr::filter(os_status != 'DECEASED')

stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months
)
rm(censored_data)
rm(observed_data)
str(stan_data)
```

(wrap this prep-data process in a function `gen_stan_data` for later)

```{r gen-stan-data-function}
gen_stan_data <- function(data) {
    observed_data <- data %>%
        dplyr::filter(os_status == 'DECEASED')
    
    censored_data <- data %>%
        dplyr::filter(os_status != 'DECEASED')
    
    stan_data <- list(
        Nobs = nrow(observed_data),
        Ncen = nrow(censored_data),
        yobs = observed_data$os_months,
        ycen = censored_data$os_months
    )
}
```

#### test simulated values with stan

Let's call `stan`:

```{r first-stan-run, warning = TRUE}
recover_simulated <- 
    rstan::stan(stan_file,
                data = gen_stan_data(simulated_data),
                chains = 4,
                iter = 1000,
                seed = 1328025050
                )
print(recover_simulated)
```


What's wrong with this picture?

 (A: poor convergence)
 (A: in some chains, we see a lot of numerical problems.)

#### Setting initial values

This step is usually optional, but may be necessary for some models. 

In this case, it may be useful to set initial values. Recall the projected range of our transformed parameter `alpha`?

By default, Stan chooses a random initial value for each parameter on the unconstrained scale between -2 and 2. This random initialization is on the *unconstrained support* for each parameter. This guarantees that initial values are consistent with the constrained range.

When we pass the initial values in, however, these are on the *constrained scale*. See the [Stan  manual](http://mc-stan.org/documentation/) for more details about transformations applied to constrained variables.

##### gen_inits function 

Let's review the parameters block for this model again.

```{r}
print_stan_file(stan_file, section = 'parameters')
```

We have two parameters for which initial values should be set. 

Let's try modifying the initial range for `alpha_raw` to utilize a smaller range than the default.

```{r stan-init-values}
gen_inits <- function() {
      list(
        alpha_raw = 0.01*rnorm(1),
        mu = rnorm(1)
      )
}
```

We wrap this in a function so that each chain will have a different set of initial values.

#### stan code with initial values

Let's try fitting our stan model again with our initial values function.

```{r sim-stanfit-with-inits}
recover_simulated2 <- 
    rstan::stan(stan_file,
                data = gen_stan_data(simulated_data),
                chains = 4,
                iter = 1000,
                init = gen_inits
                )
print(recover_simulated2)
```

Now we see fewer numerical problems, and better R-hat values.

#### checking convergence 

Assessing convergence can be a tricky business, since every model & every scenario is different. 

In general, I look at three things when assessing convergence : 

1. Rhat values (are they close to 1)?
2. Review traceplots for `lp__` & key parameters
3. launch `shinystan` for further checking.

##### Reviewing traceplots

In this case, the traceplot of `lp__` (the log-posterior) looks like it's well mixed:

```{r sim-traceplot-lp, fig.height=3}
rstan::traceplot(recover_simulated2, 'lp__')
```

Similarly, those for our parameters of interest look good:

```{r sim-traceplot-params}
rstan::traceplot(recover_simulated2, c('alpha','mu'), ncol = 1)
```

##### Launch shiny-stan

We could also launch [shinystan](http://github.com/stan-dev/shinystan) to check for divergent transitions or excessive autocorrelation in the chains. This is particularly helpful when diagnosing convergence issues. 

```{r sim-launch-shinystan, eval = F}
if (interactive())
    shinystan::launch_shinystan(recover_simulated2)
```

### Reviewing posterior distributions of parameters

Next, we can check to see whether the model was able to recover our parameters.

#### extracting parameters from the Stanfit object

We use the `rstan::extract()` function to extract parameters from the 
stanfit object.

E.g. to extract `alpha` & `mu`:

```{r sim-extract-alpha}
pp_alpha <- rstan::extract(recover_simulated2,'alpha')$alpha
pp_mu <- rstan::extract(recover_simulated2,'mu')$mu
```

Each of these is a 1xD vector of values, where D = the number of posterior (post-warmup) draws.

In this case, we have 2000: 4 chains * 1000 iterations / 2

For example, we can plot the posterior distribution of `alpha`, in context of `test_alpha`:

```{r plot-alpha-vs-test}
ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
    geom_density(aes(x = alpha)) + 
    geom_vline(aes(xintercept = test_alpha), colour = 'red') +
    ggtitle('Posterior distribution of alpha\nshowing true value in red')
```

and, repeat the same for `mu`:

```{r plot-mu-vs-test}
ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
    geom_density(aes(x = mu)) + 
    geom_vline(aes(xintercept = test_mu), colour = 'red') +
    ggtitle('Posterior distribution of mu\nshowing true value in red')
```

The posterior estimates of the parameters *do* contain those used to simulate our data, but they are not necessarily at the mode of each distribution. 

Also, we have a high degree of correlation between `mu` and `alpha`:

```{r plot-mu-vs-alpha}
ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
    geom_density2d(aes(x = alpha, y = mu)) +
    geom_point(aes(x = test_alpha, y = test_mu), colour = 'red', size = 2) +
    ggtitle('Posterior distributions of mu and alpha\nshowing true parameter values in red')
```

It does seem like these values are on the edge of the posterior density for the parameter values.

*How likely are the parameter values used to simulate our data, according to this model?*

One nice thing about sampling from the full posterior is that we can compute posterior probabilities directly, by summarizing over the values.

To compute the probability of seeing a value of `r stringr::str_c('alpha >= ',test_alpha)`:

```{r}
min(mean(pp_alpha >= test_alpha), mean(pp_alpha <= test_alpha))
```

And for `mu` (testing `r stringr::str_c('mu >= ', test_mu)`:

```{r}
min(mean(pp_mu >= test_mu), mean(pp_mu <= test_mu))
```

The joint probability of seeing this combination of parameter values, however, may not be as encouraging.

We could calculate this by, for example:

```{r}
mean(pp_mu >= test_mu & pp_alpha >= test_alpha)
```

Though, you would want to adjust the direction of `<` or `>` depending on the observed values.

### Posterior predictive checks

Next we might ask whether this error in recovering the parameters used to simulate data are substantive. Perhaps we can be a little off in estimating the baseline hazard parameters, so long as our inferences about biomarkers are sustained?

To do this, we will simulate data from our posterior draws of parameters. These are called the **posterior predicted values**. Their distribution is the **posterior predictive distribution**.

*(It may seem crazy to do this for simulated data, but it's not really extra work since we will likely want to re-use this process on our observed dataset. We will take care throughout to write each step as a function to be re-used later.)*

#### simulating data for each posterior draw

We can use hadley's `purrr::map2` to simulate data for each pair of `mu`*`alpha` values.

```{r sim-post-predict}
pp_newdata <- 
    purrr::map2(.x = pp_alpha,
                .y = pp_mu,
                .f = ~ weibull_sim_data(alpha = .x, 
                                mu = .y,
                                n = test_n
                                )
                )
```

If you're not familiar with `purrr`, I recommend you inspect the resulting object.

We now have a list of D datasets, each containing a simulation according to that draw's parameter values for `mu` & `alpha`.

Let's plot the time to event in the posterior draws, and compare this to the test dataset we used to fit our model.

```{r sim-plot-time-to-event}
ggplot(pp_newdata %>%
           dplyr::bind_rows() %>%
           dplyr::mutate(type = 'posterior predicted values') %>%
           bind_rows(simulated_data %>% dplyr::mutate(type = 'actual data'))
       , aes(x = os_months, group = os_status, colour = os_status, fill = os_status)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~type, ncol = 1)
```

These look pretty similar.

#### summarizing posterior predictive draws

Next we might ask about the posterior estimates of the survival curve. How would we estimate this?

One way (there may be several) is to:

1. compute the cumulative survival at each observed timepoint for each draw from the posterior
2. aggregate the cumulative survival estimates to discrete units of time
3. summarize the cumulative survival for each interval, over the posterior draws.

This is the method we will use here.

```{r sim-pp-survdata}
## cumulative survival rate for each posterior draw
pp_survdata <-
    pp_newdata %>%
    purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
    purrr::map(~ survival::survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
    purrr::map(fortify)

## summarize cum survival for each unit time (month), summarized at 95% confidence interval
pp_survdata_agg <- 
    pp_survdata %>%
    purrr::map(~ dplyr::mutate(., time_group = floor(time))) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(time_group) %>%
    dplyr::summarize(surv_mean = mean(surv)
                     , surv_p50 = median(surv)
                     , surv_lower = quantile(surv, probs = 0.025)
                     , surv_upper = quantile(surv, probs = 0.975)
                     ) %>%
    dplyr::ungroup()
```

Finally, we overlay the posterior predictive simulations of the survival curve with that from our original test dataset.

```{r sim-plot-ppcheck}
## km-curve for test data 
test_data_kmcurve <- 
    fortify(
        survival::survfit(
            Surv(os_months, os_deceased) ~ 1, 
            data = simulated_data %>% 
                dplyr::mutate(os_deceased = os_status == 'DECEASED')
            )) %>%
    dplyr::mutate(lower = surv, upper = surv)

ggplot(pp_survdata_agg %>%
           dplyr::mutate(type = 'posterior predicted values') %>%
           dplyr::rename(surv = surv_p50, lower = surv_lower, upper = surv_upper, time = time_group) %>%
           bind_rows(test_data_kmcurve %>% dplyr::mutate(type = 'actual data')),
       aes(x = time, group = type, linetype = type)) + 
    geom_line(aes(y = surv, colour = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    xlim(c(0, 200))
```

Here we see that the survival curve from our original (simulated) data matches the posterior predictive density pretty closely. This is a good thing since we are working with simulated data!

#### Saving as a function 

As before, we will want to wrap this in a function so that it can be reused in future steps, e.g. when we work with our TCGA data.

```{r pp_predict-function}
pp_predict_surv <- function(pp_alpha, pp_mu, n,
                            level = 0.9,
                            plot = F, data = NULL,
                            sim_data_fun = weibull_sim_data
                            ) {
    pp_newdata <- 
        purrr::map2(.x = pp_alpha,
                    .y = pp_mu,
                    .f = ~ sim_data_fun(alpha = .x, mu = .y, n = n)
                    )
    
    pp_survdata <-
        pp_newdata %>%
        purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
        purrr::map(~ survival::survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
        purrr::map(fortify)
    
    ## compute quantiles given level 
    lower_p <- 0 + ((1 - level)/2)
    upper_p <- 1 - ((1 - level)/2)
    
    pp_survdata_agg <- 
        pp_survdata %>%
        purrr::map(~ dplyr::mutate(.,
                                   time_group = floor(time))) %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(time_group) %>%
        dplyr::summarize(surv_mean = mean(surv)
                         , surv_p50 = median(surv)
                         , surv_lower = quantile(surv,
                                                 probs = lower_p)
                         , surv_upper = quantile(surv,
                                                 probs = upper_p)
                         ) %>%
        dplyr::ungroup()
    
    if (plot == FALSE) {
        return(pp_survdata_agg)
    } 
    
    ggplot_data <- pp_survdata_agg %>%
           dplyr::mutate(type = 'posterior predicted values') %>%
           dplyr::rename(surv = surv_p50,
                         lower = surv_lower,
                         upper = surv_upper, time = time_group)
    
    if (!is.null(data))
        ggplot_data <- 
            ggplot_data %>% 
            bind_rows(
                fortify(
                    survival::survfit(
                        Surv(os_months, os_deceased) ~ 1, 
                        data = data %>% 
                            dplyr::mutate(
                                os_deceased = os_status == 'DECEASED')
                        )) %>%
                dplyr::mutate(lower = surv,
                              upper = surv, type = 'actual data')
                )
    
    pl <- ggplot(ggplot_data,
                 aes(x = time, group = type, linetype = type)) + 
        geom_line(aes(y = surv, colour = type)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
        
    pl 
}
```

## Fitting our model to TCGA data

We are now ready to fit our model to the TCGA data. Here, we can reuse the several functions we defined in earlier steps.

```{r fit-model-tgca}
wei_fit <- rstan::stan(file = stan_file,
                       data = gen_stan_data(clin_data),
                       iter = 1000,
                       chains = 4,
                       init = gen_inits
                       )
print(wei_fit)
```

## Checking convergence

In the previous section, we reviewed various ways to check convergence.

### Inspect traceplots

We then also look at the traceplot of the log-posterior. I usually check this first since it gives me a sense of the overall model fit.

```{r lp-traceplot}
rstan::traceplot(wei_fit, c('lp__'), ncol = 1)
```

And, for key parameters (in this case, `alpha` & `mu`)

```{r param-traceplot}
rstan::traceplot(wei_fit, c('alpha','mu'), ncol = 1)
```

### Launch shinystan 

For more detailed model inspection we can leverage the awesome `shinystan`. 
This displays best-practice diagnostics (e.g. autocorrelation of chains) in an easy-to-use interface.

```{r launch-shinystan, eval = F}
if (interactive())
    launch_shinystan(wei_fit)
```

## Posterior predictive checks

And, we run our posterior predictive checks using the function we defined earlier:

```{r wei-ppchecks}
pl <- pp_predict_surv(pp_alpha = extract(wei_fit,'alpha')$alpha,
                pp_mu = extract(wei_fit,'mu')$mu,
                n = nrow(clin_data),
                data = clin_data,
                plot = T
                ) 
pl + 
    xlim(NA, 150) +
    ggtitle('Posterior predictive checks for NULL weibull model\nfit to TCGA data; showing 90% CI')
```

This is a pretty close fit to our data, although it's hard to say just how close of a fit it is.

*For what proportion of timepoints is the observed event rate within the 90% CI of the posterior predicted values?*

```{r summarize-coverage}
## summarize 90% CI of predicted event rate for each interval
pp_agg <- pp_predict_surv(pp_alpha = extract(wei_fit,'alpha')$alpha,
                pp_mu = extract(wei_fit,'mu')$mu,
                n = nrow(clin_data)
                )


## summarize observed data into same time_groups
act_agg <- 
    survival::survfit(Surv(os_months, I(os_status == 'DECEASED')) ~ 1,
                             data = clin_data
                             ) %>%
    fortify() %>%
    dplyr::mutate(time_group = floor(time)) %>%
    dplyr::group_by(time_group) %>%
    dplyr::summarise(observed_surv = mean(surv)) %>%
    dplyr::ungroup()

## compute proportion of observed values within 90% ci
act_agg %>%
    dplyr::inner_join(pp_agg, by = 'time_group') %>%
    dplyr::mutate(within_interval = ifelse(observed_surv >= surv_lower & observed_surv <= surv_upper,
                                           1, 0),
                  time_set = cut(time_group, breaks = c(0, 100))
                  ) %>%
    dplyr::group_by(time_set) %>%
    dplyr::summarize(mean(within_interval))
```

This isn't particularly encouraging, since it suggests that far fewer than the expected 90% of observations fall within the posterior intervals.

But, we haven't included any covariate values. Does this improve with additional covariates?

# Second analysis: Parametric survival with covariates

We have included a modified version of our NULL model code that includes 
a covariate matrix `X`. 

This model also includes a prior on the coefficient estimates which has the result of regularizing the estimates. Following the [Peltola et al, 2014](http://ceur-ws.org/Vol-1218/bmaw2014_paper_8.pdf) example, a different prior is placed on known "background" coefficients than on potential biomarkers whose correlation to the outcome of interest is less well established.

Before looking at biomarkers, however, we will first examine how well our model with the `bg` coefficients fits.

Let's review the Stan file for this model:

```{r fullmodel-stan-file}
stan_file <- system.file('stan', 'weibull_survival_model.stan', package =  'biostan')

biostan::print_stan_file(stan_file)
```

Or, if you are at a console, you might want to open the file to view its contents:

```{r fullmodel-edit-file}
## open stan file to review contents 
if (interactive())
    file.edit(stan_file)
```

### Review data block

Our data block now looks like this : 

```{r full-data-block}
biostan::print_stan_file(stan_file, section = 'data')
```

with 3 new inputs expected:

- **M_bg**: number of 'bg' covariates (background or known covariates)
- **Xobs_bg**: matrix of covariates for records with observed events
- **Xcen_bg**: matrix of covariates for records with censored events


### Review parameters block

Our parameters block now contains some new parameters, which are used to control the degree of regularization of the beta coefficients: 

```{r full-parameters-block}
biostan::print_stan_file(stan_file, section = 'parameters')
```

The regularization is done by calling a function (defined in the functions block) in `transformed parameters`:

```{r full-functions-block}
biostan::print_stan_file(stan_file, section = 'transformed parameters')
```

### Review model block

Our model block also now contains the expansion of our linear predictor, `mu`, for censored and uncensored observations.

```{r}
biostan::print_stan_file(stan_file, section = 'model')
```

## Prepare input data for stan

As before, we will need to prepare our input list to pass data to the compiled Stan model.

We will update our `gen_stan_data` function to include these covariates:

```{r fullmodel-gen-stan-data}
gen_stan_data <- function(data, formula = as.formula(~ 1)) {
    
    if (!inherits(formula, 'formula'))
        formula <- as.formula(formula)
    
    observed_data <- data %>%
        dplyr::filter(os_status == 'DECEASED')
    
    censored_data <- data %>%
        dplyr::filter(os_status != 'DECEASED')
    
    Xobs_bg <- observed_data %>%
        model.matrix(formula, data = .)
    
    Xcen_bg <- censored_data %>% 
        model.matrix(formula, data = . )
    
    assertthat::assert_that(ncol(Xcen_bg) == ncol(Xobs_bg))
    M_bg <- ncol(Xcen_bg)
    
    if (M_bg > 1) {
        if ("(Intercept)" %in% colnames(Xobs_bg))
            Xobs_bg <- array(Xobs_bg[,-1], dim = c(nrow(observed_data), M_bg - 1))
        if ("(Intercept)" %in% colnames(Xcen_bg))
            Xcen_bg <- array(Xcen_bg[,-1], dim = c(nrow(censored_data), M_bg - 1))
        assertthat::assert_that(ncol(Xcen_bg) == ncol(Xobs_bg))
        M_bg <- ncol(Xcen_bg)
    }
    
    stan_data <- list(
        Nobs = nrow(observed_data),
        Ncen = nrow(censored_data),
        yobs = observed_data$os_months,
        ycen = censored_data$os_months,
        M_bg = M_bg,
        Xcen_bg = array(Xcen_bg, dim = c(nrow(censored_data), M_bg)),
        Xobs_bg = array(Xobs_bg, dim = c(nrow(observed_data), M_bg))
    )
}
```

This function will take a `formula` object as input & construct the `stan_data` input.

For example: 

```{r fullmodel-test-gen-data}
stan_input_data <- gen_stan_data(clin_data, '~ I(gender=="MALE")')

str(stan_input_data)
```

### Testing model on sample data 

Let's try to run this on our updated stan code

```{r fullmodel-test-stan}
testfit <- rstan::stan(stan_file,
                       data = gen_stan_data(clin_data, '~ I(gender=="MALE")'),
                       init = gen_inits,
                       iter = 4,
                       chains = 1
                       )

fullfit <- rstan::stan(stan_file,
                       data = gen_stan_data(clin_data, '~ I(gender=="MALE")'),
                       init = gen_inits,
                       iter = 5000,
                       chains = 4
                       )
```

### Update inits function

We are seeing some poor sampling for certain parameters - let's try updating our inits function to reduce the range of initial values given.

We have 5 parameters of interest:

```{r fullmodel-param}
biostan::print_stan_file(stan_file, section = 'parameters')
```

This function takes a parameter (`M_bg`) and returns a function to generate initial values for each chain.

```{r gen-inits-2}
gen_inits2 <- function(M_bg) {
    function()
      list(
        alpha_raw = 0.01*rnorm(1),
        mu = rnorm(1),
        tau_s_bg_raw = 0.1*abs(rnorm(1)),
        tau_bg_raw = array(abs(rnorm(M_bg)), dim = c(M_bg)),
        beta_bg_raw = array(rnorm(M_bg), dim = c(M_bg))
      )
}
```

Let's test this out with our model. We will use a fairly innocuous covariate (`gender` as an example).

```{r fullmodel-test-with-inits}
testfit2 <- rstan::stan(stan_file,
                       data = gen_stan_data(clin_data, '~ I(gender=="MALE")'),
                       init = gen_inits2(M_bg = 1),
                       iter = 4,
                       chains = 1
                       )

fullfit2 <- rstan::stan(stan_file,
                       data = gen_stan_data(clin_data, '~ I(gender=="MALE")'),
                       init = gen_inits2(M_bg = 1),
                       iter = 2000,
                       chains = 4,
                       control = list(adapt_delta = 0.995, max_treedepth = 10)
                       )
```

#### aside on input data types

While that model is running, let's take a look back at the `init` object we created. Both this & the `gen_stan_data` functions contain the use of `array(..., dim = c(.,.))` syntax. This is due to Stan's sensitivity to the structure of inputs given to it.

For an example of how sensitive `stan()` is to dimension of inputs, try running the following exercize.

Compare the value of *a vector wrapped by `array()`*

```{r stan-alt-inits1a}
array(rnorm(0), dim = c(0))
```

with that of *a naked vector*

```{r stan-alt-inits2a}
rnorm(0)
```

These two look the same, don't they?

**However, looks can be deceiving. **

The data type stan is looking for is a `vector`: 

```{r review-stan-params}
biostan::print_stan_file(stan_file, section = 'parameters')
```

Ironically, the only way to force an R object to be a vector is to wrap it in the `array()` operator.

For example, compare the structure of *a vector wrapped in an array* 

```{r stan-alt-inits1b}
str(array(rnorm(0), dim = c(0)))
```

with the structure of *a naked vector*

```{r stan-alt-inits2b}
str(rnorm(0))
```

These are different objects in the eyes of `stan`. 

### Reviewing model convergence

Let's review the fit object from our model.

```{r fullmodel-review-output}
print(fullfit2)
```

Based on this, we can't immediately say that the model didn't converge - since our `Rhat` values are close to 1. 

We can also say that the model has a coefficient estimate for `genderMale` that is slightly less than 0. 

Let's look at the traceplots : 

```{r fullmodel-traceplot-lp}
rstan::traceplot(fullfit2, 'lp__')
```

```{r fullmodel-traceplot-alpha-mu}
rstan::traceplot(fullfit2, c('alpha','mu'), ncol = 1)
```

```{r fullmodel-traceplot-beta}
rstan::traceplot(fullfit2, 'beta_bg')
```

If we looked at shiny-stan for this model, we would see some divergent transitions (which isn't great). But these appear to have been sampled reasonably well. 

## Reviewing parameter estimates

Before spending too much time figuring out how to sample from this model's posterior distribution perfectly, let's take a look at our parameter estimates.

Our goal is to 

1. determine if the model is useful
2. check the posterior predicted values against our data

### Summarizing inference for coefficients

First let's take a look at the distribution of our `beta` estimates. What is the range of impact this parameter will likely have on outcome?

```{r fullmodel-summarize-gender}
pp_beta_bg <- rstan::extract(fullfit2, 'beta_bg')$beta_bg
ggplot(data = data.frame(beta_bg = unlist(pp_beta_bg)),
        aes(x = beta_bg)) + 
    geom_density()
```

This is not surprising, since remember our coefficient estimates are regularized using a laplace prior. 

How likely is this coefficient to be > 0 ?

```{r fullmodel-sum-beta}
mean(pp_beta_bg >= 0)
```

How well does this model fit our data ? 

```{r}
## among women
pp_predict_surv(pp_alpha = rstan::extract(fullfit2, 'alpha')$alpha,
                pp_mu = rstan::extract(fullfit2,'mu')$mu,
                n = nrow(clin_data %>% dplyr::filter(gender != 'MALE')),
                data = clin_data %>% dplyr::filter(gender != 'MALE'),
                plot = TRUE
                )

## among men
pp_predict_surv(pp_alpha = rstan::extract(fullfit2, 'alpha')$alpha,
                pp_mu = rstan::extract(fullfit2,'mu')$mu + rstan::extract(fullfit2,'beta_bg')$beta_bg,
                n = nrow(clin_data %>% dplyr::filter(gender == 'MALE')),
                data = clin_data %>% dplyr::filter(gender == 'MALE'),
                plot = TRUE
                )

```

