pseudospotBoundaries <- function(me, 
                                 distance = 100, 
                                 start_margin = 55/2, 
                                 size_of_spot = 55) {
    
    require(MoleculeExperiment)
    require(dplyr)
    require(purrr)
    
    # Updated_extent function to output the extent of each sample
    updated_extent <- function(me, 
                               assayName = "detected") {
        molucules <- molecules(me, assayName = assayName, flatten = TRUE)
        extent_per_sample <- molucules %>%
            group_by(sample_id) %>%
            summarise(xmin = min(x_location),
                      xmax = max(x_location),
                      ymin = min(y_location),
                      ymax = max(y_location))
        return(extent_per_sample)
    }
    
    # Integrate funciton for pseudo_spots
    # distance: The distance between the centers of adjacent pseudospots.
    # start_margin: The distance from the center of the first spot to the edge of the molecular map boundaries.
    # size_of_spot: The diameter of each pseudospot.
    
    # Function to create the spot based on the extent
    create_spots <- function(extent, 
                             distance = 100, 
                             start_margin = 55/2, 
                             size_of_spot = 55) {
        angles <- seq(0, 360, by = 10)
        radius <- size_of_spot / 2
        
        # Function to create spots for a single sample
        create_spots_single_sample <- function(sample_extent) {
            xmin <- sample_extent$xmin + start_margin
            xmax <- sample_extent$xmax - start_margin
            ymin <- sample_extent$ymin + start_margin
            ymax <- sample_extent$ymax - start_margin
            
            x_central_odd <- seq(from = xmin, to = xmax, by = distance)
            x_central_even <- seq(from = xmin + distance / 2, to = xmax, by = distance)
            y_central_step <- distance * sqrt(3) / 2
            y_central <- seq(from = ymax, to = ymin, by = -y_central_step)
            
            spots_odd <- expand.grid(x = x_central_odd, y = y_central[c(TRUE, FALSE)], KEEP.OUT.ATTRS = FALSE)
            spots_even <- expand.grid(x = x_central_even, y = y_central[c(FALSE, TRUE)], KEEP.OUT.ATTRS = FALSE)
            
            spots <- rbind(spots_odd, spots_even) %>%
                mutate(spot_id = row_number()) %>%
                filter(x <= xmax, x >= xmin, y <= ymax, y >= ymin)
            
            # Calculate the boundary coordinates for each spot
            spots_boundary <- map2_df(spots$x, spots$y, ~{
                tibble(
                    sample_id = sample_extent$sample_id,
                    spot_id = rep(spots$spot_id[spots$x == .x & spots$y == .y], length(angles)),
                    x_central = rep(.x, length(angles)),
                    y_central = rep(.y, length(angles)),
                    x_location = .x + radius * cos(angles * pi / 180),
                    y_location = .y + radius * sin(angles * pi / 180)
                )
            })
            
            return(spots_boundary)
        }
        
        # Apply function to each sample and combine results
        all_spots <- bind_rows(lapply(split(extent, extent$sample_id), create_spots_single_sample))
        
        return(all_spots)
    }
    
    extent = updated_extent(me, assayName = "detected")
    
    # Generate spots using the function
    spots_boundary <- create_spots(extent, distance, start_margin, size_of_spot)
    
    pseudospotMEList <- dataframeToMEList(spots_boundary,
                                          dfType = "boundaries",
                                          assayName = "pseudospot",
                                          sampleCol = "sample_id",
                                          factorCol = "spot_id",
                                          xCol = "x_location",
                                          yCol = "y_location"
    )
    
    boundaries(me, "pseudospot") <- pseudospotMEList
    
    return(me)
}
