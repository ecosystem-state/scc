scc
========
The 'scc' repository is a hub for ecosystem state indices in the Southern California Current. Files in the 
code/ folder contain the R scripts needed to replciate the analysis. Included are (1) models for just the 
climate data alone, (2) models for just the biological time series alone, and (3) models that jointly use the
climate data to predict observed biological responses. Additional details of the approach, and results, can
be found in Hunsicker et al. (2021). 

bayesdfa
========

Much of these models implements a Bayesian version of Dynamic Factor Analysis (DFA). The package 'bayesdfa'
is on CRAN, and the homepage is [https://fate-ewi.github.io/bayesdfa/articles/bayesdfa.html](https://fate-ewi.github.io/bayesdfa/articles/bayesdfa.html)  

You can install the development version of the package with:

``` r
# install.packages("devtools")
devtools::install_github("fate-ewi/bayesdfa")
```


NOAA Disclaimer
========

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

