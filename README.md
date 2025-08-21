## Introduction

This R Markdown document demonstrates how to calculate the IBS
(Identity-By-State) genetic distance between genotypes using a sample
genotype matrix. The distance calculation function is sourced from
[ibs-distance.R](ibs-distance.R).

------------------------------------------------------------------------

## Load Required Packages

We need the `reshape2` package for data transformation (specifically the
`melt()` function) to process the genetic distance matrix into a
long-form data frame.

    # Load necessary library
    library(reshape2)

------------------------------------------------------------------------

## Install IBS Distance Function

The IBS distance function `ibs.dist` is installed from this GitHub
repository:

-   **Repository**: <https://github.com/ignaciojci/ibsdiste>

<!-- -->

    # Install and load library
    if (!require("remotes")) install.packages("remotes")
    remotes::install_github("ignaciojci/ibsdist")

------------------------------------------------------------------------

## Example Genotype Matrix

For demonstration, we create a small genotype matrix with 3 genotypes
(Geno1, Geno2, Geno3) and 4 markers (Mark1 to Mark4). The values
represent the allele counts per marker.

    # Define a sample genotype matrix
    m <- matrix(c(2,0,0,
                  0,0,2,
                  0,0,2,
                  2,2,2),
                nrow=3,
                dimnames = list(paste0("Geno",1:3),
                                paste0("Mark",1:4)))

    # Display the matrix
    m

    ##       Mark1 Mark2 Mark3 Mark4
    ## Geno1     2     0     0     2
    ## Geno2     0     0     0     2
    ## Geno3     0     2     2     2
    
    # A larger matrix
    M <- matrix(sample(c(0L,1L,2L,NA_integer_), 100*500, replace=TRUE),
            nrow=2000, ncol=5000)

------------------------------------------------------------------------

## Compute IBS Genetic Distance

We compute the genetic distance matrix using the `ibs.dist()` function.

    # Calculate the IBS genetic distance matrix
    gd <- ibs.dist(m)

    # Display the resulting distance matrix
    gd

    ##       Geno1 Geno2 Geno3
    ## Geno1  0.00  0.25  0.75
    ## Geno2  0.25  0.00  0.50
    ## Geno3  0.75  0.50  0.00

------------------------------------------------------------------------

## Convert Distance Matrix to Long-Format Table

To make the pairwise distances easier to read and analyze, we reshape
the distance matrix into a long-format data frame. Only the upper
triangle of the symmetric distance matrix is used to avoid redundant
pairs.

    # Convert distance matrix into long format for easier interpretation
    gd_list <- melt(gd, varnames = c("LineA", "LineB"), value.name = "Genetic Distance")[upper.tri(gd), ]

    # Display the formatted genetic distances
    gd_list

    ##   LineA LineB Genetic Distance
    ## 4 Geno1 Geno2             0.25
    ## 7 Geno1 Geno3             0.75
    ## 8 Geno2 Geno3             0.50
