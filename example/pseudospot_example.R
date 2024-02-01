library(MoleculeExperiment)
library(ggplot2)
library(dplyr)
library(purrr)
library(ggforce) 

repoDir <- system.file("extdata", package = "MoleculeExperiment")
repoDir <- paste0(repoDir, "/xenium_V1_FF_Mouse_Brain")
me <- MoleculeExperiment::readXenium(repoDir,
                                      keepCols = "essential",
                                      addBoundaries = "cell")

source("../R/pseudospotBoundaries.R")

# identify pseudospot boundaries assuming that the x,y coordinates are in microns
 
me <- pseudospotBoundaries(me, 
                           distance = 100, 
                           start_margin = -55, 
                           size_of_spot = 55)

spot_overlay <- ggplot_me() +
    # add cell segments and colour by cell id
    geom_polygon_me(me, byFill = "segment_id", colour = "black", alpha = 0.1) +
    # add nuclei segments and colour the border with red
    geom_polygon_me(me, assayName = "pseudospot", fill = NA, colour = "red", size = 0.2) 
spot_overlay

# create a spatialExperiment object using pseudospots as boundaries
# this aggregates the gene expression within each spot
# the spatialCoords of spe now point to the centres of the pseudospots  
spe <- countMolecules(me, boundariesAssay = "pseudospot")
