library(MoleculeExperiment)
library(scater)

repoDir <- system.file("extdata", package = "MoleculeExperiment")
repoDir <- paste0(repoDir, "/xenium_V1_FF_Mouse_Brain")
me <- readXenium(repoDir,
                 keepCols = "essential",
                 addBoundaries = "cell"
)

source("../R/tileBoundaries.R")

me <- tileBoundaries(me, tile_width = 1)
extent(me)

spe = countMolecules(me, boundariesAssay = "tiles")


spe <- addPerCellQCMetrics(spe)
spe <- spe[,spe$total > 0]

plotReducedDim(spe, "spatial", colour_by = "total", other_fields = list("sample_id"),
               point_size = 5) + 
  facet_wrap(~sample_id) + 
  coord_fixed() + 
  geom_polygon(aes(x = x_location, y = y_location, group = segment_id),
               fill = NA,
               colour = "grey",
               data = MoleculeExperiment::boundaries(me,
                                 assayName = "tiles",
                                 flatten = TRUE)) +
  geom_polygon(aes(x = x_location, y = y_location, group = segment_id),
               fill = NA,
               colour = "black",
               data = MoleculeExperiment::boundaries(me,
                                                     assayName = "cell",
                                                     flatten = TRUE)) +
    geom_point(aes(x = x_location, y = y_location), data = MoleculeExperiment::molecules(me, flatten = TRUE),
             size = 0.3) + 
  NULL
