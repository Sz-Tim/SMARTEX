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
POI_bbox <- POI |> filter(name=="JC241 mooring") |>
  st_as_sf(coords=c("lon", "lat"), crs=4326) |>
  st_transform(crs="ESRI:102007") |>
  st_buffer(dist=300000) |>
  st_transform(4326) |>
  st_bbox() |>
  map(~round(.x, 2))
east_bbox <- list(xmin=-130, xmax=-90, ymin=2, ymax=20)
ccz_bbox <- list(xmin=-160, xmax=-90, ymin=0, ymax=23.5)

# bathy <- read_stars("data/bathymetry/gebco_2023_n23.5_s0.0_w-160.0_e-90.0.tif")

machine <- switch(Sys.info()["nodename"], 
                  "SA04TS-CB33BW2"=1, 
                  "SA05MH-5SJZN53"=2, 
                  "salmon"=3)
dirs <- list(
  nc=c("D:/SMARTEX/data/HYCOM/",
       # "W:/common/sa04ts/HYCOM/",
       "D:/Projects/SMARTEX/data/HYCOM/",
       "~/SMARTEX/data/HYCOM/"),
  rds=c("D:/SMARTEX/data/HYCOM/",
        "D:/Projects/SMARTEX/data/HYCOM/",
        "~/SMARTEX/data/HYCOM/"),
  rng=c("D:/SMARTEX/data/HYCOM/",
        "D:/Projects/SMARTEX/data/HYCOM/",
        "~/SMARTEX/data/HYCOM/")
) |> map(~.x[machine])

bbox <- east_bbox
date_seq <- seq(ymd("2023-07-12"), ymd("2023-07-24"), by=1)



# hycom indexes -----------------------------------------------------------

thredds_url <- "http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0"
hy_nc <- nc_open(thredds_url)
hy_i <- list(lon=ncvar_get(hy_nc, "lon"),
             lat=ncvar_get(hy_nc, "lat"),
             time=ncvar_get(hy_nc, "time"),
             depth=ncvar_get(hy_nc, "depth"))
nc_close(hy_nc)

for(j in seq_along(date_seq)) {
  library(tidyverse); library(glue)
  library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
  library(stars)
  library(sf)
  library(ncdf4)
  library(ggnewscale)
  library(colorspace)
  library(av)
  library(tictoc)
  gc()
  d <- date_seq[j]
  ymdh_j <- as.numeric(difftime(as_datetime(paste(d, "00:00:00")), 
                                as_datetime("2000-01-01 00:00:00"),
                                units="hours"))
  cat("Starting", as.character(d), "at", format(Sys.time()), "\n")
  
  file_base <- "hycom_glby_930_ccz_surface"
  
  for(i in seq(0, 23, by=3)) {
    dates_ij <- rep(ymdh_j + i, 2)
    i_ <- str_pad(i, 2, "left", "0")
    # try({
      download_hycom(thredds_url, hy_i, 
                     bbox, dates_ij, 
                     depth_range=c(1, 1), 
                     out_nc=glue("{dirs$nc}{file_base}_{d}_{i_}h.nc"),
                     out_rds=glue("{dirs$rds}{file_base}_{d}_{i_}h.rds"),
                     out_rng=glue("{dirs$rng}{file_base}_ranges_{d}_{i_}h.rds")
                     )
    # }, outFile=glue("~/HYCOM.log"))
    gc()
  }
}





# mean u, v by lon, lat, depth --------------------------------------------

# calculating mean u, v at each lon/lat/depth
cube_f <- dirf(dirs$rds, "3D.rds") |> grep("ranges", x=_, invert=T, value=T)
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
saveRDS(hycom_sum, glue("{dirs$rds}UK1_uv_means.rds"))
rm(hycom_sum); gc()


# calculating min/max/sd EKE at each depth
cube_f <- dirf(dirs$rds, "3D.rds") |> grep("ranges", x=_, invert=T, value=T)
cube_dates <- str_sub(cube_f, -21, -5)
hycom_means <- readRDS(glue("{dirs$rds}UK1_uv_means.rds"))
for(i in 1:length(cube_f)) {
  cube_i <- readRDS(cube_f[i]) |>
    select(-starts_with("EKE"), -ends_with("_mn")) |>
    left_join(hycom_means, by=c("id", "depth")) |>
    mutate(EKE=0.5*( (u-u_mn)^2 + (v-v_mn)^2) * 1e4)
  saveRDS(cube_i, cube_f[i])
  dirf(dirs$rng, glue("ranges_{cube_dates[i]}.rds")) |>
    readRDS() |>
    select(-starts_with("EKE")) |>
    left_join(cube_i |> group_by(depth) |>
                summarise(EKE_min=min(EKE, na.rm=T), 
                          EKE_max=max(EKE, na.rm=T),
                          EKE_mean=mean(EKE, na.rm=T),
                          EKE_sd=sd(EKE, na.rm=T)) |>
                ungroup()) |>
    mutate(across(everything(), 
                  ~if_else(is.infinite(.x) | is.nan(.x), NA, .x))) |> 
    saveRDS(dirf(dirs$rng, glue("ranges_{cube_dates[i]}.rds")))
  gc()
}




# plots -------------------------------------------------------------------

figs_temp <- "D:/SMARTEX/figs/temp/"
ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
hycom_means <- readRDS(glue("{dirs$rds}UK1_uv_means.rds"))
ranges <- dirf(dirs$rds, "3D.rds") |> grep("ranges", x=_, value=T) |>
  map_dfr(readRDS) |>
  group_by(depth) |>
  summarise(across(ends_with("min"), min),
            across(ends_with("max"), max),
            across(ends_with("sd"), median),
            across(ends_with("mean"), mean)) |>
  ungroup() |>
  mutate(u2v2_rng=u2v2_max - u2v2_min,
         EKE_rng=EKE_max - EKE_min) 
rng <- list(u=with(ranges, range(c(u_min, u_max), na.rm=T)),
            v=with(ranges, range(c(v_min, v_max), na.rm=T)),
            u2v2=with(ranges, range(c(u2v2_min, u2v2_max), na.rm=T)),
            uv=with(ranges, range(c(uv_min, uv_max), na.rm=T)),
            EKE=with(ranges, range(c(EKE_min, EKE_max), na.rm=T)),
            u_std=with(ranges, range(c((u_min-u_mean)/u_sd, 
                                       (u_max-u_mean)/u_sd), na.rm=T)),
            v_std=with(ranges, range(c((v_min-v_mean)/v_sd, 
                                       (v_max-v_mean)/v_sd), na.rm=T)),
            u2v2_std=with(ranges, pmin(range(c((u2v2_min-u2v2_mean)/u2v2_sd, 
                                               (u2v2_max-u2v2_mean)/u2v2_sd), na.rm=T), 10)),
            uv_std=with(ranges, pmin(range(c((uv_min-uv_mean)/uv_sd, 
                                             (uv_max-uv_mean)/uv_sd), na.rm=T), 10)),
            EKE_std=with(ranges, pmin(range(c((EKE_min-EKE_mean)/EKE_sd, 
                                              (EKE_max-EKE_mean)/EKE_sd), na.rm=T), 10)))

cube_f <- dirf(dirs$rds, "3D.rds") |> grep("ranges", x=_, invert=T, value=T)

if(FALSE) {
  readRDS(cube_f[1]) |> filter(time==first(c(time)), depth==0) |>
    select(lon, lat, id) |>
    mutate(lon=if_else(lon > 180, lon - 360, lon)) |>
    st_as_sf(coords=c("lon", "lat"), crs=4326) |>
    write_sf("data/HYCOM/HYCOM_UK1_coords.gpkg", delete_layer=TRUE) 
}
hycom_sf <- st_read("data/HYCOM/HYCOM_UK1_coords.gpkg") 
ctd_eddyTrack <- readRDS(last(dirf("data/AVISO", "tracks_ccz.*rds"))) |>
  # filter(track==132761) |>
  select(lon, lat, date, track) |>
  st_as_sf(coords=c("lon", "lat"), crs=4326) |>
  st_crop(st_bbox(hycom_sf)) |>
  arrange(track, date) |>
  group_by(track) |>
  summarise(do_union=F) |>
  st_cast("LINESTRING") |> st_segmentize(5e3) |>
  st_cast("MULTIPOINT") |> st_cast("POINT")
eddy_id <- st_nearest_feature(ctd_eddyTrack, hycom_sf) |> unique()
hycom_ids <- hycom_sf$id[eddy_id]
eddy_hycom_i <- tibble(id=hycom_ids,
                       track_id=1:length(hycom_ids))

# Parameters to plot
depth_subsets <- c(0, 500, 1000)
hycom_par <- tibble(par=c("u", "v", "u2v2", "uv", 
                          "u_std", "v_std", "u2v2_std", "uv_std",
                          "uvDir",
                          "EKE", "EKE_std"),
                    flab=c("u", "v", "u2v2", "uv",
                           "z-u", "z-v", "z-u2v2", "z-uv",
                           "uvDir", 
                           "EKE", "z-EKE"),
                    dark=c(T, T, F, F, 
                           T, T, F, F,
                           T,
                           F, F))

rm(hycom_sf); rm(bathy); rm(ctd_eddyTrack); gc()
pkg <- c("tidyverse", "glue", "terra", "sf",  "ggnewscale", "colorspace")
# library(foreach)
# library(doParallel)
# cl <- makeCluster(1)
# registerDoParallel(cl)
for(j in seq_along(cube_f)) {
  chunk_df <- map_dfr(cube_f[j], readRDS) |>
    mutate(lon=if_else(lon > 180, lon - 360, lon),
           date=as_date(time),
           hour=hour(time)) |>
    left_join(ranges |> select(depth, ends_with("sd"), ends_with("mean"))) |>
    # left_join(hycom_means, by=c("id", "depth")) |>
    mutate(EKE_std=pmin((EKE-EKE_mean)/EKE_sd, 10),
           u_std=(u-u_mean)/u_sd,
           v_std=(v-v_mean)/v_sd,
           uv_std=pmin((uv-uv_mean)/uv_sd, 10),
           u2v2_std=pmin((u2v2-u2v2_mean)/u2v2_sd, 10)) |>
    select(-ends_with("sd")) |>
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
  # tracks_east <- readRDS(last(dirf("data/AVISO", "tracks_east.*rds"))) |>
  #   filter(track==132761) |>
  #   mutate(active=factor(track==132761, levels=c("FALSE", "TRUE"))) |>
  #   filter(between(lon, lon_rng[1], lon_rng[2]))
  
  
  for(i in seq_along(timesteps)) {
  # foreach(i = seq_along(timesteps), .packages=pkg) %dopar% {
    theme_set(theme_classic())
    i_time <- timesteps[i]
    i_date <- as_date(timesteps)[i]
    i_hour <- hour(timesteps)[i]
    i_df <- chunk_df |> filter(date==i_date & hour==i_hour) 
    # i_eddy <- eddy_df |> filter(date==i_date & hour==i_hour) |>
    #   group_by(lon, depth) |> slice_head(n=1) |> ungroup()
    # i_tracks <- tracks_east |> filter(date <= i_date)
    
# depth panels ------------------------------------------------------------

    for(k in 1:nrow(hycom_par)) {
      if(hycom_par$par[k] == "uvDir") {
        plot_hycom_depthPanels(
          df=i_df, depths=depth_subsets, ccz=ccz, POI=POI,
          fill_var="uvDir", fill_lim=NA, alpha_lim=rng$uv_std,
          xlim=lon_rng, ylim=lat_rng,
          title=i_date, out_dim=c(6, 9),
          out_f=glue("{figs_temp}UK1_uvDir_",
                     "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
          tracks=NULL, darkLines=T)
      } else {
        plot_hycom_depthPanels(
          df=i_df, depths=depth_subsets, ccz=ccz, POI=POI,
          fill_var=hycom_par$par[k], fill_lim=rng[[hycom_par$par[k]]],
          xlim=lon_rng, ylim=lat_rng,
          title=i_date, out_dim=c(6, 9),
          out_f=glue("{figs_temp}UK1_{hycom_par$flab[k]}_",
                     "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
          tracks=NULL, darkLines=hycom_par$dark[k])
      }
      gc()
    }
    
    
    
# eddy path ---------------------------------------------------------------

    # for(d in c(500, 4000)) {
    #   for(k in 1:nrow(hycom_par)) {
    #     if(hycom_par$par[k] == "uvDir") {
    #       plot_hycom_eddyProfile(
    #         df=i_eddy, depthLim=d, lon_diff=lon_diff, POI=POI[4,],
    #         fill_var="uvDir", fill_lim=NA, alpha_lim=rng$uv_std,
    #         xlim=lon_rng, ylim=lat_rng,
    #         title=i_date, out_dim=c(9, 7),
    #         out_f=glue("{figs_temp}eddy{d}m_uvDir_",
    #                    "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
    #         tracks=i_tracks, darkLines=T)
    #     } else {
    #       plot_hycom_eddyProfile(
    #         df=i_eddy, depthLim=d, lon_diff=lon_diff, POI=POI[4,],
    #         fill_var=hycom_par$par[k], fill_lim=rng[[hycom_par$par[k]]],
    #         xlim=lon_rng, ylim=lat_rng,
    #         title=i_date, out_dim=c(9, 7),
    #         out_f=glue("{figs_temp}eddy{d}m_{hycom_par$flab[k]}_",
    #                    "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
    #         tracks=i_tracks, darkLines=hycom_par$dark[k])
    #     }
    #     gc()
    #   }
    # }
  }
}
# stopCluster(cl)




# animate -----------------------------------------------------------------

library(tictoc); library(av); library(sevcheck); library(glue)
img_path <- "figs/anim/temp"
img_path <- "D:/SMARTEX/figs/temp/"
out_path <- "figs/anim_temp/"

fps <- 8*5
sets <- c(outer(c("UK1_"),#, "eddy4000m_", "eddy2000m_", "eddy1000m_", "eddy500m_"), 
                c("u", "v", "u2v2", "uv", 
                  "z-u", "z-v", "z-u2v2", "z-uv",
                  "uvDir", "EKE", "z-EKE"), 
                paste0))

for(i in sets) {
  try({
    cat("starting", i, "\n")
    tic()
    av::av_encode_video(dirf(img_path, paste0(i, "_")), 
                        glue("{out_path}{i}.mp4"),
                        framerate=fps)
    toc()
  }, silent=T)
  gc()
}








poi_df <- dirf("D:/Projects/SMARTEX/data/HYCOM", "POI_cube_2023-0[3-5].*rds") |>
  map_dfr(~readRDS(.x) |> 
            filter(abs(c(lon)-243.52) < 0.01, 
                   abs(c(lat)-13.89) < 0.03))

poi_df |> 
  filter(depth < 4500) |> 
  ggplot(aes(ymin=depth0, ymax=depth1, xmin=time-1.5*60*60, xmax=time+1.5*60*60, fill=EKE)) + 
  geom_rect() + 
  scale_fill_viridis_c("HYCOM\nEKE", option="turbo", begin=0.1) + 
  scale_y_reverse() + 
  labs(x="Date", y="Depth (m)", title="JC241 mooring location") + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"))
ggsave("figs/JC241_EKE_2023-MarMay_HYCOM.png", width=10, height=4)

poi_df |> 
  filter(depth <= 500) |> 
  ggplot(aes(ymin=depth0, ymax=depth1, xmin=time-1.5*60*60, xmax=time+1.5*60*60, fill=EKE)) + 
  geom_rect() + 
  scale_fill_viridis_c("HYCOM\nEKE", option="turbo", begin=0.1) + 
  scale_y_reverse() + 
  labs(x="Date", y="Depth (m)", title="JC241 mooring location") + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"))
ggsave("figs/JC241_EKE_0-500m_2023-MarMay_HYCOM.png", width=10, height=4)

poi_df |> 
  filter(depth < 4500) |> 
  group_by(depth) |>
  mutate(EKE=c(scale(EKE))) |> 
  ungroup() |> 
  ggplot(aes(ymin=depth0, ymax=depth1, xmin=time-1.5*60*60, xmax=time+1.5*60*60, fill=EKE)) + 
  geom_rect() + 
  scale_fill_viridis_c("HYCOM\nz(EKE)", option="turbo", begin=0.1) + 
  scale_y_reverse() + 
  labs(x="Date", y="Depth (m)", title="JC241 mooring location") + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"))
ggsave("figs/JC241_EKE-z_2023-MarMay_HYCOM.png", width=10, height=4)

poi_df |> 
  filter(depth <= 500) |> 
  group_by(depth) |>
  mutate(EKE=c(scale(EKE))) |> 
  ungroup() |> 
  ggplot(aes(ymin=depth0, ymax=depth1, xmin=time-1.5*60*60, xmax=time+1.5*60*60, fill=EKE)) + 
  geom_rect() + 
  scale_fill_viridis_c("HYCOM\nz(EKE)", option="turbo", begin=0.1) + 
  scale_y_reverse() + 
  labs(x="Date", y="Depth (m)", title="JC241 mooring location") + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"))
ggsave("figs/JC241_EKE-z_0-500m_2023-MarMay_HYCOM.png", width=10, height=4)




poi_df |> 
  filter(depth < 4500) |> 
  ggplot(aes(ymin=depth0, ymax=depth1, xmin=time-1.5*60*60, xmax=time+1.5*60*60, fill=uv)) + 
  geom_rect() + 
  scale_fill_viridis_c("HYCOM\nuv", option="turbo", begin=0.1) + 
  scale_y_reverse() + 
  labs(x="Date", y="Depth (m)", title="JC241 mooring location") + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"))
ggsave("figs/JC241_uv_2023-MarMay_HYCOM.png", width=10, height=4)

poi_df |> 
  filter(depth <= 500) |> 
  ggplot(aes(ymin=depth0, ymax=depth1, xmin=time-1.5*60*60, xmax=time+1.5*60*60, fill=uv)) + 
  geom_rect() + 
  scale_fill_viridis_c("HYCOM\nuv", option="turbo", begin=0.1) + 
  scale_y_reverse() + 
  labs(x="Date", y="Depth (m)", title="JC241 mooring location") + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"))
ggsave("figs/JC241_uv_0-500m_2023-MarMay_HYCOM.png", width=10, height=4)

poi_df |> 
  filter(depth < 4500) |> 
  group_by(depth) |>
  mutate(uv=c(scale(uv))) |> 
  ungroup() |> 
  ggplot(aes(ymin=depth0, ymax=depth1, xmin=time-1.5*60*60, xmax=time+1.5*60*60, fill=uv)) + 
  geom_rect() + 
  scale_fill_viridis_c("HYCOM\nz(uv)", option="turbo", begin=0.1) + 
  scale_y_reverse() + 
  labs(x="Date", y="Depth (m)", title="JC241 mooring location") + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"))
ggsave("figs/JC241_uv-z_2023-MarMay_HYCOM.png", width=10, height=4)

poi_df |> 
  filter(depth <= 500) |> 
  group_by(depth) |>
  mutate(uv=c(scale(uv))) |> 
  ungroup() |> 
  ggplot(aes(ymin=depth0, ymax=depth1, xmin=time-1.5*60*60, xmax=time+1.5*60*60, fill=uv)) + 
  geom_rect() + 
  scale_fill_viridis_c("HYCOM\nz(uv)", option="turbo", begin=0.1) + 
  scale_y_reverse() + 
  labs(x="Date", y="Depth (m)", title="JC241 mooring location") + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"))
ggsave("figs/JC241_uv-z_0-500m_2023-MarMay_HYCOM.png", width=10, height=4)
