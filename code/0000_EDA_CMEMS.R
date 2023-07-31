


rm(list=ls()); gc()

library(terra)
library(tidyverse)
library(lubridate)
library(sf)
library(ncdf4)
library(oce)
library(gganimate)
library(ggnewscale)
theme_set(theme_classic())



ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
bathy_reclass <- cbind(c(-Inf, seq(-6000, 0, by=500)),
                       c(seq(-6000, 0, by=500), Inf),
                       seq(-6.5, 0, by=0.5)*1e3)
bathy <- rast("data/bathymetry/gebco_2023_n23.5_s0.0_w-160.0_e-110.0.tif") |>
  aggregate(fact=20) |>
  classify(bathy_reclass) |>
  as.polygons() |> 
  st_as_sf() |> 
  rename_with(~"elev", .cols=1)
ccz_bbox <- list(xmin=360-160, xmax=360-111.4067, ymin=0, ymax=23.5)
cmems_nc <- nc_open("data/CMEMS/dataset-duacs-nrt-global-merged-allsat-phy-l4_1689690749440.nc")

cmems_ls <- list(time=c(ncvar_get(cmems_nc, "time")),
                 lon=c(ncvar_get(cmems_nc, "longitude")),
                 lat=c(ncvar_get(cmems_nc, "latitude")),
                 adt=ncvar_get(cmems_nc, "adt"),
                 sla=ncvar_get(cmems_nc, "sla"))
nc_close(cmems_nc)

cmems_df <- expand_grid(time=cmems_ls$time,
                        lat=cmems_ls$lat,
                        lon=cmems_ls$lon) |>
  mutate(time=ymd("1950-01-01") + time,
         adt=c(cmems_ls$adt),
         sla=c(cmems_ls$sla)) |>
  mutate(sla_disc=case_when(sla < -0.05 ~ -1,
                            sla > -0.05 & sla < 0.05 ~ 0,
                            sla > 0.05 ~ 1))

anim <- cmems_df |>
  # filter(time < "2022-06-08") |>
  filter(sla_disc != 0) |>
  ggplot() + 
  geom_sf(data=bathy, aes(fill=elev), colour=NA) +
  scale_fill_distiller(type="seq", palette="Greys") +
  new_scale_fill() +
  geom_tile(aes(lon, lat, fill=sla), alpha=0.5) + 
  scale_fill_distiller("sla", type="div", palette="RdBu") +
  geom_sf(data=ccz, linewidth=0.25, colour="grey30") +
  transition_manual(time) + 
  labs(title="{closest_state}") + 
  theme(axis.title=element_blank())
anim_save("figs/anim/00_sla_test.gif", anim, fps=20,
          width=8, height=4, units="in", res=300)
