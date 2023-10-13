# HYCOM EDA
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

POI <- read_csv("data/SMARTEX_locations_copy.csv") |>
  filter(grepl("BGR|Long|Atoll|001|002", Name_Full))

# ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
POI_bbox <- list(xmin=-120, xmax=-100, ymin=7, ymax=15)
east_bbox <- list(xmin=-130, xmax=-90, ymin=2, ymax=20)
ccz_bbox <- list(xmin=-160, xmax=-90, ymin=0, ymax=23.5)

# bathy <- read_stars("data/bathymetry/gebco_2023_n23.5_s0.0_w-160.0_e-90.0.tif")

machine <- switch(Sys.info()["nodename"], 
                  "SA04TS-CB33BW2"=1, 
                  "SA05MH-5SJZN53"=2, 
                  "salmon"=3)
dirs <- list(
  nc=c(#"E:/Projects/SMARTEX/SMARTEX/data/HYCOM/",
       "W:/common/sa04ts/HYCOM/",
       "D:/Projects/SMARTEX/data/HYCOM/",
       "~/SMARTEX/data/HYCOM/"),
  rds=c("E:/Projects/SMARTEX/SMARTEX/data/HYCOM/",
        "D:/Projects/SMARTEX/data/HYCOM/",
        "~/SMARTEX/data/HYCOM/"),
  rng=c("E:/Projects/SMARTEX/SMARTEX/data/HYCOM/",
        "D:/Projects/SMARTEX/data/HYCOM/",
        "~/SMARTEX/data/HYCOM/")
) |> map(~.x[machine])


hy_dataset <- c("uv3z", "ssh")[2]
bbox <- ccz_bbox
date_seq <- seq(ymd("2023-01-01"), today(), by=1)



# hycom indexes -----------------------------------------------------------

thredds_url <- "http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0"
hy_nc <- nc_open(thredds_url)
hy_i <- list(lon=ncvar_get(hy_nc, "lon"),
             lat=ncvar_get(hy_nc, "lat"),
             time=ncvar_get(hy_nc, "time"),
             depth=ncvar_get(hy_nc, "depth"))
nc_close(hy_nc)

for(j in seq_along(date_seq)) {
  gc()
  d <- date_seq[j]
  ymdh_j <- as.numeric(difftime(as_datetime(paste(d, "00:00:00")), 
                                as_datetime("2000-01-01 00:00:00"),
                                units="hours"))
  cat("Starting", as.character(d), "at", format(Sys.time()), "\n")
  
  file_base <- glue("hycom_glby_930_{hy_dataset}_",
                    "{bbox$xmin}Eto{bbox$xmax}E_", 
                    "{bbox$ymin}Nto{bbox$ymin}N")
  
  for(i in seq(0, 23, by=3)) {
    dates_ij <- rep(ymdh_j + i, 2)
    i_ <- str_pad(i, 2, "left", "0")
    if(file.exists(glue("{dirs$nc}{file_base}_{d}_{i_}h.nc"))) {
      next
    }
    try({
      download_hycom(glue("{thredds_url}/{hy_dataset}"), hy_i, 
                     bbox, dates_ij, 
                     depth_range=c(1, ifelse(hy_dataset=="ssh", 1, 40)), 
                     out_nc=glue("{dirs$nc}{file_base}_{d}_{i_}h.nc"),
                     out_rds=glue("{dirs$rds}{file_base}_{d}_{i_}h.rds"),
                     out_rng=glue("{dirs$rng}{file_base}_ranges_{d}_{i_}h.rds"))
    }, outFile=glue("~/HYCOM.log"))
  }
  gc()
}





# mean u, v by lon, lat, depth --------------------------------------------

# calculating mean u, v at each lon/lat/depth
cube_f <- dirf(glue("{E_dir}data/HYCOM"), "POI_cube_2023")
hycom_sum <- readRDS(cube_f[1]) |>
    group_by(id, depth) |> 
    summarise(u_sum=sum(u, na.rm=T), 
              v_sum=sum(v, na.rm=T),
              u_n=sum(!is.na(u)),
              v_n=sum(!is.na(v))) |> 
    ungroup()
for(i in 2:length(cube_f)) {
  cube_i <- readRDS(cube_f[i]) |>
    group_by(id, depth) |> 
    summarise(u_sum=sum(u, na.rm=T), 
              v_sum=sum(v, na.rm=T),
              u_n=sum(!is.na(u)),
              v_n=sum(!is.na(v))) |> 
    ungroup()
  hycom_sum[,3:6] <- hycom_sum[,3:6] + cube_i[,3:6]
}
hycom_sum <- hycom_sum |>
  mutate(u_mn=u_sum/u_n,
         v_mn=v_sum/v_n) |>
  select(-ends_with("_sum"), -ends_with("_n"))
saveRDS(hycom_sum, glue("{E_dir}data/HYCOM/POI_uv_means.rds"))
rm(hycom_sum); gc()


# calculating min/max/sd EKE at each depth
cube_f <- dirf(glue("{E_dir}data/HYCOM"), "POI_cube_2023")
cube_dates <- str_sub(cube_f, -14, -5)
hycom_means <- readRDS(glue("{E_dir}data/HYCOM/POI_uv_means.rds"))
for(i in 1:length(cube_f)) {
  cube_i <- readRDS(cube_f[i]) |>
    select(-starts_with("EKE"), -ends_with("_mn")) |>
    left_join(hycom_means, by=c("id", "depth")) |>
    mutate(EKE=0.5*( (u-u_mn)^2 + (v-v_mn)^2) * 1e4)
  saveRDS(cube_i, cube_f[i])
  readRDS(glue("{E_dir}data/HYCOM/temp/POI_ranges_{cube_dates[i]}.rds")) |>
    select(-starts_with("EKE")) |>
    left_join(cube_i |> group_by(depth) |>
                summarise(EKE_min=min(EKE, na.rm=T), 
                          EKE_max=max(EKE, na.rm=T),
                          EKE_mean=mean(EKE, na.rm=T),
                          EKE_sd=sd(EKE, na.rm=T)) |>
                ungroup()) |>
    saveRDS(glue("{E_dir}data/HYCOM/temp/POI_ranges_{cube_dates[i]}.rds"))
  gc()
}




# plots -------------------------------------------------------------------

hycom_means <- readRDS(glue("{E_dir}data/HYCOM/POI_uv_means.rds"))
ranges <- dirf(glue("{E_dir}data/HYCOM/temp"), "POI_ranges_") |>
  map_dfr(readRDS) |>
  group_by(depth) |>
  summarise(across(ends_with("min"), min),
            across(ends_with("max"), max),
            across(ends_with("sd"), median)) |>
  ungroup() |>
  mutate(u2v2_rng=u2v2_max - u2v2_min)
rng <- list(u=with(ranges, range(c(u_min, u_max))),
            v=with(ranges, range(c(v_min, v_max))),
            u2v2=with(ranges, range(c(u2v2_min, u2v2_max))),
            uv=with(ranges, range(c(uv_min, uv_max))),
            EKE=with(ranges, range(c(EKE_min, EKE_max))),
            u_std=with(ranges, range(c(u_min/u_sd, u_max/u_sd))),
            v_std=with(ranges, range(c(v_min/v_sd, v_max/v_sd))),
            u2v2_std=c(0,1),
            uv_std=with(ranges, range(c(uv_min/uv_sd, uv_max/uv_sd))),
            EKE_std=with(ranges, range(c(EKE_min/EKE_sd, EKE_max/EKE_sd))))

cube_f <- dirf(glue("{E_dir}data/HYCOM"), "POI_cube_2023")

if(FALSE) {
  readRDS(cube_f[1]) |> filter(time==first(time), depth==0) |>
    select(lon, lat, id) |>
    st_as_sf(coords=c("lon", "lat"), crs=4326) |>
    write_sf("data/HYCOM/HYCOM_POI_coords.gpkg", delete_layer=TRUE) 
}
hycom_sf <- st_read("data/HYCOM/HYCOM_POI_coords.gpkg") 
ctd_eddyTrack <- readRDS(last(dirf("data/AVISO", "tracks_ccz.*rds"))) |>
  filter(track==132761) |>
  select(lon, lat, date) |>
  st_as_sf(coords=c("lon", "lat"), crs=4326) |>
  st_crop(st_bbox(hycom_sf)) |>
  arrange(date) |>
  summarise(do_union=F) |>
  st_cast("LINESTRING") |> st_segmentize(5e3) |>
  st_cast("MULTIPOINT") |> st_cast("POINT")
eddy_id <- st_nearest_feature(ctd_eddyTrack, hycom_sf) |> unique()
hycom_ids <- hycom_sf$id[eddy_id]
eddy_hycom_i <- tibble(id=hycom_ids,
                       track_id=1:length(hycom_ids))

# Parameters to plot
depth_subsets <- c(0, 500, 1000, 1500, 2000)
hycom_par <- tibble(par=c(#"u", "v", "u2v2", "uv", 
                          #"u_std", "v_std", "u2v2_std", "uv_std",
                          "uvDir",
                          "EKE", "EKE_std"),
                    flab=c(#"u", "v", "u2v2", "uv",
                           #"z-u", "z-v", "z-u2v2", "z-uv",
                           "uvDir", 
                           "EKE", "z-EKE"))

rm(hycom_sf); rm(bathy); rm(ctd_eddyTrack); gc()
pkg <- c("tidyverse", "glue", "terra", "sf",  "ggnewscale", "colorspace")
library(foreach)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
for(j in seq_along(cube_f)) {
  chunk_df <- map_dfr(cube_f[j], readRDS) |>
    mutate(lon=if_else(lon > 180, lon - 360, lon),
           date=as_date(time),
           hour=hour(time)) |>
    left_join(ranges |> select(depth, ends_with("sd"), starts_with("u2v2"))) |>
    left_join(hycom_means, by=c("id", "depth")) |>
    mutate(EKE_std=EKE/EKE_sd,
           u_std=u/u_sd,
           v_std=v/v_sd,
           uv_std=uv/uv_sd,
           u2v2_std=(u2v2-u2v2_min)/u2v2_rng) |>
    select(-ends_with("sd"), -u2v2_min, -u2v2_max, -u2v2_rng) |>
    ungroup()
  eddy_df <- chunk_df |>
    right_join(eddy_hycom_i)
  if(!"track_id" %in% names(POI)) {
    POI <- POI %>%
      mutate(track_id=eddy_df$track_id[st_nearest_feature(
        st_as_sf(POI, coords=c("lon", "lat"), crs=4326),
        st_as_sf(eddy_df, coords=c("lon", "lat"), crs=4326))]
      )
  }
  lon_diff <- mean(diff(sort(unique(chunk_df$lon))))
  timesteps <- sort(unique(chunk_df$time))
  lon_rng <- range(chunk_df$lon)
  lat_rng <- range(chunk_df$lat)
  tracks_east <- readRDS(last(dirf("data/AVISO", "tracks_east.*rds"))) |>
    filter(track==132761) |>
    mutate(active=factor(track==132761, levels=c("FALSE", "TRUE"))) |>
    filter(between(lon, lon_rng[1], lon_rng[2]))
  
  
  
  foreach(i = seq_along(timesteps), .packages=pkg) %dopar% {
    theme_set(theme_classic())
    i_time <- timesteps[i]
    i_date <- as_date(timesteps)[i]
    i_hour <- hour(timesteps)[i]
    i_df <- chunk_df |> filter(date==i_date & hour==i_hour) 
    i_eddy <- eddy_df |> filter(date==i_date & hour==i_hour) |>
      group_by(lon, depth) |> slice_head(n=1) |> ungroup()
    i_tracks <- tracks_east |> filter(date <= i_date)
    
# depth panels ------------------------------------------------------------

    for(k in 1:nrow(hycom_par)) {
      if(hycom_par$par[k] == "uvDir") {
        plot_hycom_depthPanels(
          df=i_df, depths=depth_subsets, ccz=ccz, POI=POI,
          fill_var="uvDir", fill_lim=NA, alpha_lim=rng$uv_std,
          xlim=lon_rng, ylim=lat_rng,
          title=i_date, out_dim=c(6, 9),
          out_f=glue("figs/anim/temp/POI_uvDir_",
                     "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
          tracks=i_tracks)
      } else {
        plot_hycom_depthPanels(
          df=i_df, depths=depth_subsets, ccz=ccz, POI=POI,
          fill_var=hycom_par$par[k], fill_lim=rng[[hycom_par$par[k]]],
          xlim=lon_rng, ylim=lat_rng,
          title=i_date, out_dim=c(6, 9),
          out_f=glue("figs/anim/temp/POI_{hycom_par$flab[k]}_",
                     "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
          tracks=i_tracks)
      }
      gc()
    }
    
    
    
# eddy path ---------------------------------------------------------------

    for(d in c(1000, 4000)) {
      for(k in 1:nrow(hycom_par)) {
        if(hycom_par$par[k] == "uvDir") {
          plot_hycom_eddyProfile(
            df=i_eddy, depthLim=d, lon_diff=lon_diff, POI=POI[2,],
            fill_var="uvDir", fill_lim=NA, alpha_lim=rng$uv_std,
            xlim=lon_rng, ylim=lat_rng,
            title=i_date, out_dim=c(9, 7),
            out_f=glue("figs/anim/temp/eddy{str_sub(d,1,1)}k_uvDir_",
                       "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
            tracks=i_tracks)
        } else {
          plot_hycom_eddyProfile(
            df=i_eddy, depthLim=d, lon_diff=lon_diff, POI=POI[2,],
            fill_var=hycom_par$par[k], fill_lim=rng[[hycom_par$par[k]]], 
            xlim=lon_rng, ylim=lat_rng,
            title=i_date, out_dim=c(9, 7),
            out_f=glue("figs/anim/temp/eddy{str_sub(d,1,1)}k_{hycom_par$flab[k]}_",
                       "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
            tracks=i_tracks)
        }
        gc()
      }
    }
  }
}
stopCluster(cl)




# animate -----------------------------------------------------------------

library(tictoc); library(av); library(sevcheck); library(glue)
img_path <- "figs/anim/temp"
out_path <- "figs/anim/"

fps <- 8*5
sets <- c(outer(c("POI_", "eddy4k_", "eddy2k_", "eddy1k_"), 
                c("u", "v", "u2v2", "uv", 
                  "z-u", "z-v", "z-u2v2", "z-uv",
                  "uvDir", "EKE", "z-EKE"), 
                paste0))

for(i in sets) {
  try({
    cat("starting", i, "\n")
    tic()
    av::av_encode_video(dirf(img_path, paste0(i, "_")), 
                        glue("figs/anim/{i}.mp4"),
                        framerate=fps)
    toc()
  }, silent=T)
  gc()
}





