# Sea surface height
# SMARTEX
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk



# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(stars)
library(sf)
library(ncdf4)
library(ggnewscale)
library(colorspace)
library(av)
library(tictoc)

theme_set(theme_bw() + theme(panel.grid=element_blank()))

dirf("code/fn/") |> walk(source)
reload <- TRUE
POI <- read_csv("data/SMARTEX_locations_copy.csv") |>
  filter(grepl("BGR|Long|Atoll|001|002", Name_Full))

ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
east_bbox <- list(xmin=-130, xmax=-90, ymin=2, ymax=20)
ccz_bbox <- list(xmin=-160, xmax=-90, ymin=0, ymax=23.5)
bathy <- read_stars("data/bathymetry/gebco_2023_n23.5_s0.0_w-160.0_e-90.0.tif")
land <- st_read("data/CCZ_areas/CCZ_box_coastline.gpkg")


# Images: Eastern focus ---------------------------------------------------

if(reload) {
  cmems_east <- download_cmems(creds=readRDS("data/cmems_cred.rds"),
                               east_bbox) |>
    mutate(EKE=pmin(EKE_cm2s2, 3e3))
saveRDS(cmems_east, glue("data/CMEMS/cmems_east_{min(cmems_east$date)}_{max(cmems_east$date)}.rds"))
} else {
  cmems_east <- readRDS(last(dirf("data/CMEMS", "cmems_east")))
}


# # EKE 2D maps
# EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
#                 out_f=c("mean", "sd", "variance"),
#                 legend=c("EKE mean", "EKE sd", "EKE variance"))
# for(k in 1:nrow(EKE_i)) {
#   plot_EKE_2D(cmems_east, bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
#               east_bbox[1:2], east_bbox[3:4], 
#               glue("figs/anim/east_EKE_{EKE_i$out_f[k]}.png"), c(10, 5.5))
# }
# EKE_sf <- cmems_east |> 
#   group_by(lat, lon) |> 
#   summarise(EKE_mn=mean(EKE), 
#             EKE_var=var(EKE),
#             EKE_sd=sd(EKE)) |> 
#   ungroup() |>
#   st_as_sf(coords=c("lon", "lat"), crs=4326) 
# walk(EKE_i$var, ~EKE_sf |> select(all_of(.x), geometry) |>
#        stars::st_rasterize() |>
#        stars::write_stars(glue("data/CMEMS/{.x}_east.tif")))
#   
timesteps <- sort(unique(cmems_east$date))
cmems_vars <- list("sla"=range(cmems_east$sla, na.rm=T),
                   "EKE"=range(cmems_east$EKE, na.rm=T),
                   # "uv"=range(cmems_east$uv, na.rm=T))
                   "uv"=c(0,1))
tracks_east <- readRDS(last(dirf("data/AVISO", "tracks_east.*rds"))) |>
  filter(type=="anticyclonic")

for(i in seq_along(timesteps)) {
  i_time <- timesteps[i]
  i_df <- cmems_east |> filter(date==i_time)
  i_tracks <- tracks_east |> filter(date <= i_time) |>
    group_by(track) |>
    mutate(active=factor(track==132761, levels=c("FALSE", "TRUE"))) |>
    ungroup()

  for(k in seq_along(cmems_vars)) {
    plot_cmems(df=i_df |> filter(!is.na(sla)), bathy=bathy, land=land, ccz=ccz, POI=POI,
               fill_var=names(cmems_vars)[k], fill_lim=cmems_vars[[k]],
               alpha_lim=range(sqrt(abs(cmems_east$sla)), na.rm=T),
               xlim=east_bbox[1:2], ylim=east_bbox[3:4],
               title=i_time, out_dim=c(10, 6),
               out_f=glue("figs/anim/temp/east_{names(cmems_vars)[k]}_{i_time}.png"),
               tracks=i_tracks, darkLines=names(cmems_vars)[k]=="sla")
    gc()
  }
}
rm(cmems_east); gc()




# Images: Full CCZ --------------------------------------------------------

# cmems_ccz <- download_cmems(creds=readRDS("data/cmems_cred.rds"), ccz_bbox) |>
#   mutate(EKE=pmin(EKE_cm2s2, 3e3))
# 
# 
# # EKE 2D maps
# EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
#                 out_f=c("mean", "sd", "variance"),
#                 legend=c("EKE mean", "EKE sd", "EKE variance"))
# for(k in 1:nrow(EKE_i)) {
#   plot_EKE_2D(cmems_ccz, bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
#               ccz_bbox[1:2], ccz_bbox[3:4], 
#               glue("figs/anim/FullCCZ_EKE_{EKE_i$out_f[k]}.png"), c(11, 5))
# }
# 
# 
# timesteps <- sort(unique(cmems_ccz$date))
# sla_rng <- range(cmems_ccz$sla, na.rm=T)
# EKE_rng <- range(cmems_ccz$EKE, na.rm=T)
# tracks_ccz <- readRDS(last(dirf("data/AVISO", "tracks_ccz.*rds"))) |>
#   filter(type=="anticyclonic")
# 
# for(i in seq_along(timesteps)) {
#   i_time <- timesteps[i]
#   i_df <- cmems_ccz |> filter(date==i_time)
#   i_tracks <- tracks_ccz |> filter(date <= i_time) |>
#     group_by(track) |>
#     # mutate(active=factor(max(date)==i_time, levels=c("FALSE", "TRUE"))) |>
#     mutate(active=factor(track==132761, levels=c("FALSE", "TRUE"))) |>
#     ungroup()
# 
#   plot_cmems(df=i_df, bathy=bathy, ccz=ccz, POI=POI,
#              fill_var="sla", fill_lim=sla_rng, title=i_time,
#              out_dim=c(10, 4),
#              out_f=glue("figs/anim/temp/CCZ_sla_{i_time}.png"),
#              tracks=i_tracks)
#   gc()
#   plot_cmems(df=i_df, bathy=bathy, ccz=ccz, POI=POI,
#              fill_var="EKE", fill_lim=EKE_rng, title=i_time,
#              out_dim=c(10, 4),
#              out_f=glue("figs/anim/temp/CCZ_EKE_{i_time}.png"),
#              tracks=i_tracks)
#   gc()
# 
# }
# rm(cmems_ccz); gc()



# Images: Eastern focus, multi-year ---------------------------------------

# cmems_east <- download_cmems(creds=readRDS("data/cmems_cred.rds"), east_bbox,
#                              dates=c(ymd("2019-12-01"), today())) |>
#   mutate(EKE=pmin(EKE_cm2s2, 3e3))
# 
# # EKE 2D maps
# EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
#                 out_f=c("mean", "sd", "variance"),
#                 legend=c("EKE mean", "EKE sd", "EKE variance"))
# for(k in 1:nrow(EKE_i)) {
#   plot_EKE_2D(cmems_ccz, bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
#               ccz_bbox[1:2], ccz_bbox[3:4], 
#               glue("figs/anim/FullCCZ_EKE_{EKE_i$out_f[k]}.png"), c(11, 5))
# }
# 
# timesteps <- sort(unique(cmems_east$date))
# sla_rng <- range(cmems_east$sla, na.rm=T)
# EKE_rng <- range(cmems_east$EKE, na.rm=T)
# 
# for(i in seq_along(timesteps)) {
#   i_time <- timesteps[i]
#   i_df <- cmems_east |> filter(date==i_time)
# 
#   plot_cmems(df=i_df, bathy=bathy, ccz=ccz, POI=POI,
#              fill_var="sla", fill_lim=sla_rng,
#              xlim=east_bbox[1:2], ylim=east_bbox[3:4],
#              title=i_time, out_dim=c(10, 4),
#              out_f=glue("figs/anim/temp/eastFull_sla_{i_time}.png"))
#   gc()
#   plot_cmems(df=i_df, bathy=bathy, ccz=ccz, POI=POI,
#              fill_var="EKE", fill_lim=EKE_rng,
#              xlim=east_bbox[1:2], ylim=east_bbox[3:4],
#              title=i_time, out_dim=c(10, 4),
#              out_f=glue("figs/anim/temp/eastFull_EKE_{i_time}.png"))
#   gc()
# 
# }
# rm(cmems_east); gc()




# Images: Full CCZ, multi-year --------------------------------------------

if(reload) {
  cmems_ccz <- download_cmems(creds=readRDS("data/cmems_cred.rds"), ccz_bbox,
                              dates=c(ymd("2019-12-01"), today())) |>
    mutate(EKE=pmin(EKE_cm2s2, 3e3))
  saveRDS(cmems_ccz, glue("data/CMEMS/cmems_ccz_{min(cmems_east$date)}_{max(cmems_east$date)}.rds"))  
} else {
  cmems_ccz <- readRDS(first(dirf("data/CMEMS", "cmems_ccz")))
}

# EKE 2D maps
# EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
#                 out_f=c("mean", "sd", "variance"),
#                 legend=c("EKE mean", "EKE sd", "EKE variance"))
# for(k in 1:nrow(EKE_i)) {
#   plot_EKE_2D(cmems_ccz, bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
#               ccz_bbox[1:2], ccz_bbox[3:4], 
#               glue("figs/anim/FullCCZ_EKE_{EKE_i$out_f[k]}.png"), c(11, 5))
# }
# EKE_sf <- cmems_ccz |> 
#   group_by(lat, lon) |> 
#   summarise(EKE_mn=mean(EKE), 
#             EKE_var=var(EKE),
#             EKE_sd=sd(EKE)) |> 
#   ungroup() |>
#   st_as_sf(coords=c("lon", "lat"), crs=4326) 
# walk(EKE_i$var, ~EKE_sf |> select(all_of(.x), geometry) |>
#        stars::st_rasterize() |>
#        stars::write_stars(glue("data/CMEMS/{.x}_CCZ.tif")))

# Animations
timesteps <- sort(unique(cmems_ccz$date))
sla_rng <- range(cmems_ccz$sla, na.rm=T)
EKE_rng <- range(cmems_ccz$EKE, na.rm=T)
uv_rng <- c(0,1)

for(i in seq_along(timesteps)) {
  i_time <- timesteps[i]
  i_df <- cmems_ccz |> filter(date==i_time)
  
  plot_cmems(df=i_df |> filter(!is.na(sla)), bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="sla", fill_lim=sla_rng,
             alpha_lim=range(sqrt(abs(cmems_ccz$sla)), na.rm=T),
             title=i_time, out_dim=c(10, 4.75),
             out_f=glue("figs/anim/temp/CCZFull_sla_{i_time}.png"), darkLines=T)
  gc()
  # plot_cmems(df=i_df, bathy=NULL, ccz=ccz, POI=POI,
  #            fill_var="EKE", fill_lim=EKE_rng,  alpha_lim=c(0, max(abs(sla_rng))),
  #            title=i_time,
  #            out_dim=c(10, 4.75),
  #            out_f=glue("figs/anim/temp/CCZFull_EKE_{i_time}.png"), darkLines=F)
  # gc()
  # plot_cmems(df=i_df, bathy=NULL, ccz=ccz, POI=POI,
  #            fill_var="uv", fill_lim=uv_rng,  alpha_lim=c(0, max(abs(sla_rng))),
  #            title=i_time,
  #            out_dim=c(10, 4.75),
  #            out_f=glue("figs/anim/temp/CCZFull_uv_{i_time}.png"), darkLines=F)
  # gc()

}
rm(cmems_ccz); gc()



# Animate -----------------------------------------------------------------

img_path <- "figs/anim/temp"
out_path <- "figs/anim/"

fps <- 25
sets <- c(outer(c("CCZ_", "east_", "CCZFull_", "eastFull_"), 
                c("sla", "EKE", "uv"), 
                paste0))

for(i in sets) {
  tic()
  try({
    av_encode_video(dirf(img_path, i), 
                    glue("figs/anim/{i}.mp4"), 
                    framerate=fps)
  })
  toc()
  gc()
}
