
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1035939.svg)](https://doi.org/**.****/zenodo.********)

[![minimal R version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kotzeb0912)](https://cran.r-project.org/package=kotzeb0912) [![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-orange.svg?style=flat-square)](commits/master)

[![Last-changedate](https://img.shields.io/badge/last%20change-2021--01--27-yellowgreen.svg)](/commits/master)

# POP
## An R package in support of publication, "Species Density Models from Opportunistic Citizen Science Data." 

#### Jay M. Ver Hoef<sup>a</sup>, Devin Johnson<sup>a</sup>, Robyn Angliss<sup>a</sup>, Matt Higham<sup>b</sup>

#### <sup>a</sup>NOAA Fisheries (NMFS) Alaska Fisheries Science Center, Marine Mammal Laboratory, Seattle, WA
#### <sup>b</sup>Department of Statistics, St. Lawrence University, Canton, NY

As a scientific work, and in keeping with common scientific practicies, I kindly request that you cite my research project and applicable publications if you use my work(s) or data in your publications or presentations. Additionally, I strongly encourage and welcome collaboration to promote use of these data in the proper context and scope.  The publication is currently submitted:

#### Ver Hoef, J.M., Johnson, D., Angliss, R., and Higham, M. 2021. Species Density Models from Opportunistic Citizen Science Data. In submission.


Abstract
-----------------

 * With the advent of technology for data-gathering and storage, opportunistic citizen-science data are proliferating. Species distribution models (SDMs) aim to use species occurrence or abundance for ecological insights, prediction, and management. We analyzed a massive opportunistic data set with over 100,000 records of incidental shipboard observations of marine mammals. Our overall goal was to create maps of species density from massive opportunistic data by using spatial regression for count data with an effort offset. We illustrate the method with two marine mammals in the Gulf of Alaska and Bering Sea.

* We counted the total number of animals in 11,424 hexagons based on presence-only data. To decrease bias, we first estimated a spatial density surface for ship-days, which was our proxy variable for effort. We used spatial considerations to create pseudo-absences, and left some hexagons as missing values. Next, we created SDMs that used modeled effort to create pseudo-absences, and included the effort surface as an offset in a second stage analysis of two example species, northern fur seals and Steller sea lions.

* For both effort and species counts, we used spatial negative-binomial regression with random effects that had a multivariate normal distribution with a conditional autoregressive (CAR) covariance matrix, providing 2 million Markov chain Monte Carlo (MCMC) samples (2000 were retained) from the posterior distribution. We used a novel MCMC scheme that main tained sparse precision matrices for observed and missing data when batch sampling from the multivariate normal distribution. We also used a truncated normal distribution to stabilize estimates, and used a look-up table for sampling the autocorrelation parameter. These innovations allowed us to draw several million samples in just a few hours.

* From the posterior distributions of the SDMs, we computed two functions of interest. We normalized the SDMs and then applied an overall abundance estimate obtained from the literature to derive spatially explicit abundance estimates, especially within subsetted areas. We also created “certain hotspots” that scaled local abundance by standard deviation and using thresholds. Hexagons with values above a threshold were deemed as hotspots with enough evidence to be certain about them.

Installation
------------

Installation of this R data package is done through the `devtools::install_github()` function or by downloading the [source package from the latest release](https://github.com/jayverhoef/POP).

```
library("devtools")
install_github("jayverhoef/POP")
```

Note that some of the files in the package are over 100 Mb because they store Markov Chain Monte Carlo (MCMC) results for 11,000+ spatial polygons for 2,000 MCMC iterations.  It may take a while to install.

Data
-------------

All of the raw data used in the manuscript are contained in `POPhexagons`, which is an `sp` object of class `SpatialPolygonsDataFrame`.  After installing the package, the data can be accessed by 

```
library(POP)
data(POPhexagons)
```

To extract the `data.frame` containing the data, 

```
dataDF = POPhexagons@data
```

The `data.frame` has 5 columns, one with the polygon ID, the x- and y-coordinates of the centroids of the polygons, and the counts for shipdays, northern fur seals (labelled CU), and Steller sea lions (labelled EJ).  Missing data are denoted by `NA`.

Run R Scripts
-------------

There are 4 R scripts that were used for the data analysis.  The first is `preliminaries.R`.  This creates a `list` of each hexagon's neighbors, it creates the grid of $\rho$ values, the sparse matrices for the inverse of the covariance matrix for each $\rho$ value, and the determinants for each row table; all of which are used as a lookup table when fitting models with Markov Chain Monte Carlo (MCMC).  These are already computed and stored as data objects, so it is not necessary to run this again, but you can see how these objects were created.

```
system.file("scripts/preliminaries.R", package = "POP")
```

which created these objects stored as data:

```
data(rhogrid)
data(SigiList)
data(detsig)
data(nb_list)
data(W)
data(Dmat)
```

To run MCMC on the data, first fit the model for effort.  The script is found here

```
system.file("scripts/MCMC_shipdays.R", package = "POP")
```

This uses the function `MCMC`, and the documentation is found by

```
help(MCMC)
```
There are a few data manipulations prior to running MCMC, and these are commented.  It takes several hours for `MCMC_shipdays.R` to run, so all results have been stored.  The script itself saves results in a `list`, but this was simplified to have several data objects, which can be accessed with

```
data(MCMC_eff_beta0)
data(MCMC_eff_miss)
data(MCMC_eff_samp)
data(MCMC_eff_rho)
data(MCMC_eff_sigCAR)
data(MCMC_eff_accept_beta0)
data(MCMC_eff_accept_miss)
data(MCMC_eff_accept_samp)
data(MCMC_eff_accept_rho)
data(MCMC_eff_accept_sigCAR)
```

See documentation on `MCMC` function for fuller description.

After fitting the effort model, an MCMC model for each of the species can be run.  The scripts is found here

```
system.file("scripts/MCMC_CU.R", package = "POP")
system.file("scripts/MCMC_EJ.R", package = "POP")
```

These fits use the same function `MCMC`.

```
help(MCMC)
```
There are a few data manipulations prior to running MCMC, and these are commented in the scripts.  It takes several hours for these scripts to run, so all results have been stored.  The script itself saves results in a `list`, but this was simplified to have several data objects, which, for northern fur seals, can be accessed with

```
data(MCMC_CU_beta0)
data(MCMC_CU_miss)
data(MCMC_CU_samp)
data(MCMC_CU_rho)
data(MCMC_CU_sigma)
data(MCMC_CU_sigCAR)
data(MCMC_CU_accept_beta0)
data(MCMC_CU_accept_miss)
data(MCMC_CU_accept_samp)
data(MCMC_CU_accept_rho)
data(MCMC_CU_accept_sigma)
data(MCMC_CU_accept_sigCAR)
```

and for Steller sea lions can be accessed with 

```
data(MCMC_EJ_beta0)
data(MCMC_EJ_miss)
data(MCMC_EJ_samp)
data(MCMC_EJ_rho)
data(MCMC_EJ_sigma)
data(MCMC_EJ_sigCAR)
data(MCMC_EJ_accept_beta0)
data(MCMC_EJ_accept_miss)
data(MCMC_EJ_accept_samp)
data(MCMC_EJ_accept_rho)
data(MCMC_EJ_accept_sigma)
data(MCMC_EJ_accept_sigCAR)
```

Create Graphics
-------------

All graphs in the manuscript can be re-created using a script found here,

```
system.file("doc/figures/POPmanu_Figures.R", package = "POP")
```

Latex Document
-------------

A folder containing the latex document, and all necessary files and subfolders, that were used to create the manuscript, can be found here

```
system.file("doc", package = "POP")
```

-------------
##### Disclaimer

<sub>This repository is a scientific product and is not official communication of the Alaska Fisheries Science Center, the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All AFSC Marine Mammal Laboratory (AFSC-MML) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. AFSC-MML has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.</sub>

