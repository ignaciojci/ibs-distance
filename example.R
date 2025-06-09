# Load necessary libraries
library(reshape2)

# Source the ibs-distance.R script from GitHub
source("https://raw.githubusercontent.com/ignaciojci/ibs-distance/refs/heads/main/ibs-distance.R")

# Create the genotype matrix
m <- matrix(c(2,0,0,0,0,2,0,0,2,2,2,2), nrow=3,dimnames = list(paste0("Geno",1:3), paste0("Mark",1:4)))

# Display the matrix
m

# Calculate the IBS genetic distance
gd <- ibs.dist(m, center=TRUE)

# Melt the distance matrix into a long format and keep only the upper triangle
gd_list <- melt(gd, varnames = c("LineA", "LineB"), value.name = "Genetic Distance")[upper.tri(gd),]

# Display the resulting list
gd_list