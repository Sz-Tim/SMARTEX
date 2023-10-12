


rm(list=ls()); gc()

library(terra)
library(tidyverse)
library(sf)
library(ncdf4)
library(oce)
library(gganimate)
library(ggnewscale)
library(colorspace)
library(glue)
theme_set(theme_classic())

walk(dirf("code/fn"), source)

ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
ccz_bbox <- st_bbox(ccz) |> as.list()
bathy_reclass <- cbind(c(-Inf, seq(-6000, 0, by=500)),
                       c(seq(-6000, 0, by=500), Inf),
                       seq(-6.5, 0, by=0.5)*1e3)
bathy <- rast("data/bathymetry/gebco_2023_n23.5_s0.0_w-160.0_e-90.0.tif") |>
  aggregate(fact=20) |>
  classify(bathy_reclass) |>
  as.polygons() |> 
  st_as_sf() |> 
  rename_with(~"elev", .cols=1)

# Global Ocean Gridded L 4 Sea Surface Heights And Derived Variables Nrt
# SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046
# https://doi.org/10.48670/moi-00149
cmems_nc <- nc_open("data/CMEMS/dataset-duacs-nrt-global-merged-allsat-phy-l4_1695210696043.nc")

cmems_ls <- list(time=c(ncvar_get(cmems_nc, "time")),
                 lon=c(ncvar_get(cmems_nc, "longitude")),
                 lat=c(ncvar_get(cmems_nc, "latitude")),
                 sla=ncvar_get(cmems_nc, "sla"))
nc_close(cmems_nc)

cmems_df <- expand_grid(time=cmems_ls$time,
                        lat=cmems_ls$lat,
                        lon=cmems_ls$lon) |>
  mutate(time=ymd("1950-01-01") + time,
         sla=c(cmems_ls$sla)) |>
  mutate(sla_disc=case_when(sla < -0.05 ~ -1,
                            sla > -0.05 & sla < 0.05 ~ 0,
                            sla > 0.05 ~ 1))

rm(cmems_nc); rm(cmems_ls); gc()

cmems_temp <- cmems_df |>
  filter(sla_disc != 0)
  # filter(time < "2020-04-10")
  # filter(time < "2021-01-04")

rm(cmems_df); gc()

timesteps <- sort(unique(cmems_temp$time))
sla_rng <- range(cmems_temp$sla)
for(i in seq_along(timesteps)) {
  i_pad <- str_pad(i, 4, "left", "0")
  i_time <- timesteps[i]
  p <- cmems_temp |>
    filter(time == i_time) |>
    ggplot() + 
    geom_sf(data=bathy, aes(fill=elev), colour=NA, alpha=0.5) +
    scale_fill_distiller(type="seq", palette="Greys", direction=1) +
    new_scale_fill() +
    geom_tile(aes(lon, lat, fill=sla), alpha=0.5) + 
    # scale_fill_gradient2("sla", low="#b2182b", high="#2166ac", limits=sla_rng) +
    # scale_fill_distiller(type="div", palette="RdBu") +
    scale_fill_continuous_divergingx(palette='RdBu', mid=0, limits=sla_rng, 
                                     l3 = 0, p3 = .8, p4 = .6, rev=TRUE) + 
    geom_sf(data=ccz, linewidth=0.25, colour="grey30") +
    labs(title=i_time) + 
    theme(axis.title=element_blank())
  ggsave(glue("figs/anim/temp/CMEMS_sla-{i_pad}.png"), 
         p, width=10, height=4, units="in", dpi=200)
  gc()
}


rm(cmems_temp); gc()

# source("code/0000_make_hycom_gifs.R")


# anim <- cmems_temp |>
#   ggplot() + 
#   geom_sf(data=bathy, aes(fill=elev), colour=NA, alpha=0.5) +
#   scale_fill_distiller(type="seq", palette="Greys") +
#   new_scale_fill() +
#   geom_tile(aes(lon, lat, fill=sla), alpha=0.5) + 
#   # scale_fill_distiller("sla", type="div", palette="RdBu") +
#   scale_fill_gradient2("sla", low="#b2182b", high="#2166ac") +
#   geom_sf(data=ccz, linewidth=0.25, colour="grey30") +
#   transition_manual(time) + 
#   labs(title="{current_frame}") + 
#   theme(axis.title=element_blank())
# anim_save("figs/anim/00_sla_full.gif", anim, 
#           nframes=n_distinct(cmems_temp$time), fps=50,
#           width=10, height=4, units="in", res=150)
