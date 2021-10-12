# McCOILR

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build status](https://github.com/OJWatson/McCOILR/workflows/R-CMD-check/badge.svg)](https://github.com/OJWatson/McCOILR/actions)
[![Codecov test coverage](https://codecov.io/gh/OJWatson/McCOILR/branch/master/graph/badge.svg)](https://codecov.io/gh/OJWatson/McCOILR?branch=master)
[![Documentation via pkgdown](https://github.com/OJWatson/McCOILR/raw/master/tools/pkgdownshield.png)](https://ojwatson.github.io/McCOILR/)

### What is this?

*McCOILR* is an Rcpp wrapper for [THE REAL McCOIL](https://github.com/Greenhouse-Lab/THEREALMcCOIL) software 
developed by the [Greenhouse Lab](https://github.com/Greenhouse-Lab) that estimates complexity of infection 
and population allele frequencies using SNP data obtained from Sequenom or similar types of SNP assays. It was
simply created to aid incorporating COI estimation more easily within distributed computing pipelines, and I
claim no ownership over the original source code, and all attribution and acknowledgement should be referred to
the original [project](https://github.com/Greenhouse-Lab/THEREALMcCOIL), and the [associated publication](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005348). I have simply provided
this wrapper in the hope that it might be helpful for others, and have provided a tutorial for basic use.

***
> To view the tutorial please click [here](https://ojwatson.github.io/McCOILR/articles/introduction.html).

***

### Installing *McCOILR*

To install the development version from github the package [*devtools*](https://github.com/hadley/devtools) is required.

In order to install devtools you need to make sure you have a working development environment:

1. **Windows**: Install **[Rtools](https://cran.r-project.org/bin/windows/Rtools/)**. For help on how to install **Rtools** please see the following [guide](https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows), paying particular attention to the section about adding Rtools to your system `PATH`. 

In order to find out which version of **Rtools** you will need to check which version of R you are running. This can be be found out using the `sessionInfo()` function:

``` r 
> sessionInfo()
R version 3.3.1 (2016-06-21)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1
```

2. **Mac**: Install Xcode from the Mac App Store.

3. **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).

Once a working development environment is ready, then devtools can be installed from CRAN. If you have not set up which CRAN mirror to use
before, then you will be asked to choose your CRAN mirror. Select the cloud mirror.

```r
install.packages("devtools")
library(devtools)
```
Once devtools is installed it is best to restart our R session. To do this either close RStudio or restart R (ctrl + shift + F10). Once your R session
has been restarted the package can be installed and loaded using:

```r
devtools::install_github("OJWatson/McCOILR")
library(McCOILR)
```

***


#### Asking a question

For bug reports, feature requests, contributions, use github's [issue system.](https://github.com/OJWatson/McCOILR/issues)
