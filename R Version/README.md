TemporalNetworkAccessibility
================

-   [Overview](#overview)
-   [Package installation and Code download](#package-installation-and-code-download)
-   [Example temporal network in R](#example-temporal-network-in-r)
-   [Plot the cumulative and shortest paths](#plot-the-cumulative-and-shortest-paths)


# Overview

This software is an *R code* for the analysis of temporal networks.

It provides code used for computing the results of the Papers

**Unfolding Accessibility Provides a Macroscopic Approach to Temporal Networks**,  
Lentz et al., Phys. Rev. Lett. 110, 118701 (2013)

and

**Disease Spread through Animal Movements: A Static and Temporal Network Analysis of Pig Trade in Germany**,
Lentz et al., PLOS ONE 11, e0155196–32 (2016)

Please cite the second reference (PLOS ONE) if you use the R code.


The current version is limited to the calculation of the cumulative path duration distribution and shortest path duration distrubution in a temporal network. It consists of 4 functions:

- ```prepare_data()```: A function to prepare he data in order to set up the Sparse matices
- ```create_temporal_network()```: A function to set up a list with Sparse matrices
- ```cumul_path()```: A function to calculate the cumulative and shortest paths
- ```plot_tn()```: A function to plot the cumulative and shortest paths

A full description and the *Python* code can be found at [TemporalNetworkAccessibility](https://github.com/hartmutlentz/TemporalNetworkAccessibility)



# Package installation and Code download

## Required packages
To run the functions you need to install the following libraries 

``` r
install.packages("Matrix", dependencies = TRUE)
install.packages("pracma", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
``` 

Then load the packages ```ggplot2``` ,``` pracma``` and ```Matrix``` 

``` r
library(Matrix)
library(pracma)
library(ggplot2)
``` 

## Code download

Download the file *Unfold_accessibility.R* from github

In your script, use the functions by adding

``` r
source("Unfold_accessibility.R")
```


# Example temporal network in R

## load packages and Code

``` r
library(Matrix)
library(pracma)
library(ggplot2)
source("Unfold_accessibility.R")
``` 

## read data
Read the example dataset.

``` r
df <- read.delim("https://raw.githubusercontent.com/hartmutlentz/TemporalNetworkAccessibility/master/edgelists/sociopatterns_hypertext.dat", header = FALSE) 
```

## prepare data
Prepare the dataset.

``` r
prep_df <- prepare_dataframe(df)
```

## Create Sparse matrix
Create a list of sparse matrices (One matrix per time step)

``` r
the_file <- create_temporal_network(prep_df)
```

## Calculate cumulative paths and shortest paths durations

``` r
cp <- cumul_path(the_file)
``` 

# Plot cumulative paths and shortest paths durations 

``` r
plot_tn(cp)
``` 
