#!/usr/bin/Rscript

## Script name: bioregions.R
##
## Purpose of script:
## Cretan arthropods bioregionalisation
## using the infomap algorithm
##
## How to run:
##
## Authors: Giannis Bolanakis, Leonidas Maroulis, Savvas Paragkamian
##
## Date Created: 2024-07-17
##

library(tidyverse)
library(sf)
library(terra)
library(quadtree)
library(vegan)
library(units)


crete_shp <- sf::st_read("../data/crete/crete.shp")


#locations_inland <- data

### 1km over raster
crete_manual_grid <- rast(ncol=311,
          nrow=152,
          xmin=5523000,
          xmax=5832000,
          ymin=1400000,
          ymax=1552000,
          crs=crs(grid_1km_o))

values(crete_manual_grid) <- 0
### WGS 84 of 1km grid EEA
crete_manual_grid_w <- terra::project(crete_manual_grid, "WGS84")
values(crete_manual_grid_w) <- 0

grid_occ_count_w <- rasterize(locations_inland, crete_manual_grid_w, fun='count')
grid_occ_count_w_df <- terra::as.data.frame(grid_occ_count_w, xy=TRUE, cells=TRUE)  

######################### Adaptive resolution with quadtree ######################
qt_sampling <- quadtree(grid_sampling_count_y,
               min_cell_length=1000,
               max_cell_length=10000,
               split_threshold=4,
               split_if_any_na=F,
               split_method = "range")
qt <- qt_sampling
#qt_species <- quadtree(grid_occ_count_y,
#               min_cell_length=2000,
#               max_cell_length=10000,
#               split_threshold=20,
#               split_if_any_na=F,
#               split_method = "range")

#qt_rast <- as_raster(qt)
#qt_rast_w <- project(qt_rast, crs(locations_inland))
#qt_rast_df <- terra::as.data.frame(qt_rast_w, xy=TRUE, cells=TRUE) 

# transform the QT object directly to dataframe to avoid errors from 
# the multiple transformation from terra to sf
qt_df <- as_data_frame(qt)

# function to create polygon from extend
f1 <- function(x) {st_polygon(list(rbind(c(x[['xmin']],x[['ymin']]),
                                          c(x[['xmax']],x[['ymin']]),
                                          c(x[['xmax']],x[['ymax']]),
                                          c(x[['xmin']],x[['ymax']]),
                                          c(x[['xmin']],x[['ymin']]))))
}

qt_df$geom <- apply(X=qt_df,1 , FUN=f1)

# create the sf polygon object from the dataframe
qt_sf <- qt_df |>
    mutate(geometry=st_sfc(geom)) |>
    st_as_sf() |>
    dplyr::select(-c(geom)) |>
    st_set_crs(st_crs(locations_eea))

qt_sf_sp_base <- qt_sf |>
    st_transform(crs(locations_inland)) |> 
    st_join(locations_inland, left=F) |>
    distinct(scientificName,geometry,id)

qt_sf_sp <- qt_sf |>
    st_transform(crs(locations_inland)) |> 
    st_join(locations_inland, left=F) |>
    distinct(scientificName,geometry,id) |>
    group_by(geometry,id) |>
    mutate(n_species=n()) |>
    left_join(endemic_species, by=c("scientificName"="scientificName"))


### quadstrees plots
png(file="../plots/quadtree_crete.png",
    width = 40,
    height = 20,
    res=300,
    units = "cm",
    bg="white")
plot(qt,
     border_lwd = .01,
     xlim = c(5523000, 5834000),
     ylim = c(1420000,1548000),
     main = "expand raster dimensions")
plot(vect(grid_10km),
     lwd = .05,
     add=T)
plot(crete_eea, add=T, color="")
#plot(locations_eea_sampling, add=T, color="#d55e00")
dev.off()

### quads and biodiversity
qt_sf_all <- qt_sf_sp |>
    group_by(geometry, id) |>
    summarise(n_species=n(), .groups="keep") |>
    mutate(Order="All")

print("the summary of quads per size-area")
table(ceiling(st_area(qt_sf_all)/1000000))

qt_sf_order <- qt_sf_sp |> 
    group_by(geometry, order) |> 
    summarise(n_species=n(), .groups="keep") |>
    ungroup()

qt_sf_n_order <- qt_sf_sp |> 
    distinct(geometry, order) |> 
    group_by(geometry) |> 
    summarise(n_orders=n())

crete_quads_endemics_map <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(locations_grid, mapping=aes()) +
    geom_sf(qt_sf_all, mapping=aes(fill=n_species),
            alpha=0.8,
            colour="transparent",
            na.rm = F,
            show.legend=T) +
    scale_fill_gradient(low="#f0e442",
                        high="#d55e00",
                        guide = "colourbar")+
    coord_sf(crs="wgs84") +
    guides(fill = guide_colourbar(ticks = F,
                                  label = T,
                                  title="# endemics",
                                  title.vjust = 0.8,
                                  order = 1))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../plots/crete_quads_endemics_map.png", 
       plot=crete_quads_endemics_map, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

crete_quads_endemics_orders_map <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(qt_sf_order, mapping=aes(fill=n_species),
            alpha=0.8,
            colour="transparent",
            na.rm = FALSE,
            show.legend=T) +
    scale_fill_gradient(low="#F0E442",
                        high="#D55E00",
                        guide = "colourbar")+
    coord_sf(crs="WGS84") +
    guides(fill = guide_colourbar(ticks = FALSE,
                                  label = TRUE,
                                  title="# endemics",
                                  title.vjust = 0.8,
                                  order = 1))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank()) +
    facet_wrap(vars(order), ncol=2, scales = "fixed")

ggsave("../plots/crete_quads_endemics_orders_map.png",
       plot=crete_quads_endemics_orders_map, 
       height = 50, 
       width = 35,
       dpi = 300, 
       units="cm",
       device="png")

crete_quads_orders_map <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_sf(qt_sf_n_order, mapping=aes(fill=n_orders),
            alpha=0.8,
            colour="transparent",
            na.rm = F,
            show.legend=T) +
    scale_fill_gradient(low="#f0e442",
                        high="#d55e00",
                        guide = "colourbar")+
    coord_sf(crs="wgs84") +
    guides(fill = guide_colourbar(ticks = F,
                                  label = T,
                                  title="# orders",
                                  title.vjust = 0.8,
                                  order = 1))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../plots/crete_quads_orders_map.png", 
       plot=crete_quads_orders_map, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

