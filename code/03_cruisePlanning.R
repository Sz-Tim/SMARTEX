# Cruise planning 2024 Feb
# SMARTEX
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(ncdf4)
library(sf)
library(terra)
library(gstat)
library(colorspace)
library(av)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
dirf("code/fn/") |> walk(source)

forecast <- F
plot_altimetry <- T
hull <- F

dir_figTemp <- "D:/SMARTEX/figs/temp/"
dir_gpkg <- "data/cruise2024/"
east_bbox <- list(xmin=-120, xmax=-90, ymin=7, ymax=20)

POI <- read_csv("data/SMARTEX_locations_copy.csv") |>
  filter(grepl("BGR|Long|Atoll|001|002|006", Name_Full))
ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg") 
ccz_east <- st_read("data/CCZ_areas/ccz_outline.gpkg") |>
  st_crop(c(xmin=-120, xmax=-90, ymin=7, ymax=20))
eez <- st_read("data/CCZ_areas/CCZ_EEZs.shp") 
bathy_reclass <- cbind(c(-Inf, seq(-6000, 0, by=500)),
                       c(seq(-6000, 0, by=500), Inf),
                       seq(-6.5, 0, by=0.5)*1e3)
bathy <- rast("data/bathymetry/gebco_2023_n23.5_s0.0_w-160.0_e-90.0.tif") |>
  aggregate(fact=20) |>
  classify(bathy_reclass) |>
  as.polygons() |> 
  st_as_sf() |> 
  rename_with(~"elev", .cols=1) |>
  st_crop(unlist(east_bbox))
land <- st_read("data/CCZ_areas/CCZ_box_coastline.gpkg")



# Get CMEMS ---------------------------------------------------------------

# Global Ocean Physics Analysis and Forecast (1/12 degree)
if(forecast) {
  # https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/download?dataset=cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m
  cur_nc <- dirf("data/cruise2024", "cur_anfc.*nc") |> nc_open()
  # https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/download?dataset=cmems_mod_glo_phy_anfc_0.083deg_P1D-m
  sla_nc <- dirf("data/cruise2024", "phy_anfc.*nc") |> nc_open()
  cmems_east <- expand_grid(time=c(ncvar_get(cur_nc, "time")),
                            lat=c(ncvar_get(cur_nc, "latitude")),
                            lon=c(ncvar_get(cur_nc, "longitude"))) |>
    mutate(date=as_date(ymd("1970-01-01") + dseconds(time)),
           sla=c(ncvar_get(sla_nc, "zos")),
           uo=c(ncvar_get(cur_nc, "uo")),
           vo=c(ncvar_get(cur_nc, "vo")),
           uv=sqrt(uo^2 + vo^2),
           uvDir=atan2(vo, uo)) |>
    select(lon, lat, date, sla, uv, uvDir) |>
    arrange(date, lat, lon)
  nc_close(cur_nc); nc_close(sla_nc) 
} else {
  cmems_nc <- dirf("data/cruise2024", "obs-sl.*nc") |> nc_open()
  cmems_east <- expand_grid(time=c(ncvar_get(cmems_nc, "time")),
                            lat=c(ncvar_get(cmems_nc, "latitude")),
                            lon=c(ncvar_get(cmems_nc, "longitude"))) |>
    mutate(date=as_date(ymd("1970-01-01") + dseconds(time)),
           sla=c(ncvar_get(cmems_nc, "sla")),
           u=c(ncvar_get(cmems_nc, "ugos")),
           v=c(ncvar_get(cmems_nc, "vgos")),
           uv=sqrt(u^2 + v^2),
           uvDir=atan2(v, u),
           EKE=pmin(0.5*(c(ncvar_get(cmems_nc, "ugosa")*100)^2 + 
                           c(ncvar_get(cmems_nc, "vgosa")*100)^2),
                    2000)) |>
    select(lon, lat, date, sla, EKE, uv, uvDir) |>
    arrange(date, lat, lon)
  nc_close(cmems_nc)
}




# Get Tracks ---------------------------------------------------------------

if(hull) {
  eddy_sf <- st_read("data/cruise2024/sla-q95-hull.gpkg")
} else {
  eddy_sf <- st_read("data/cruise2024/cruise2024_eddies_SSH_peak.gpkg")
}

timesteps <- sort(unique(cmems_east$date))
timesteps <- timesteps[length(timesteps)+(-6:0)]
cmems_vars <- list("sla"=range(cmems_east$sla, na.rm=T),
                   "uv"=range(cmems_east$uv, na.rm=T))
if(!forecast) {
  cmems_vars$EKE <- range(cmems_east$EKE, na.rm=T)
}

cmems_gpkg_layers <- st_layers(glue("{dir_gpkg}cmems_{ifelse(forecast, 'Fcst', 'Obs')}.gpkg"))
for(i in seq_along(timesteps)) {
  i_time <- timesteps[i]
  i_df <- cmems_east |> 
    filter(date==i_time) |>
    mutate(id=row_number())
  i_eddy <- eddy_sf |>
    filter(date <= i_time) |>
    mutate(active=date==i_time)
  
  if(i==1) {
    cmems_grid <- i_df |>
      select(lon, lat) |>
      st_as_sf(coords=c("lon", "lat"), crs=4326) |>
      st_make_grid(n=c(n_distinct(i_df$lon), n_distinct(i_df$lat)),
                   offset=c(min(i_df$lon), min(i_df$lat)) - (1/8)) |>
      st_as_sf() |>
      mutate(id=row_number())
    cmems_grid_hiRes <- expand_grid(
      x=seq(min(i_df$lon), max(i_df$lon), by=1/16),
      y=seq(min(i_df$lat), max(i_df$lat), by=1/16)
    )
    cmems_rast <- rast(cmems_grid_hiRes, type="xyz")
  }
  
  fitsla <- gstat::gstat(formula = sla ~ 1,
                         id="sla",
                         locations=~x+y,
                         data = i_df |> rename(x=lon, y=lat) |> filter(!is.na(sla)),
                         nmax = 7,
                         set = list(idp = 2))
  if(!file.exists(glue("{dir_gpkg}interp/cmems_{ifelse(forecast, 'Fcst', 'Obs')}_sla-interp_{i_time}.tif"))) {
    cmems_rast_hiRes <- interpolate(cmems_rast, model=fitsla, index=1)
    writeRaster(cmems_rast_hiRes, 
                glue("{dir_gpkg}interp/cmems_{ifelse(forecast, 'Fcst', 'Obs')}_sla-interp_{i_time}.tif"),
                overwrite=T)
  }
  if(! i_time %in% cmems_gpkg_layers$name) {
    full_join(cmems_grid, i_df |> select(-lon, -lat), by="id") |>
      write_sf(glue("{dir_gpkg}cmems_{ifelse(forecast, 'Fcst', 'Obs')}.gpkg"), layer=i_time, append=F) 
  }
  
  if(plot_altimetry) {
    for(k in seq_along(cmems_vars)) {
      plot_cmems(df=i_df |> filter(!is.na(sla)), bathy=NULL, land=land, ccz=ccz, eez=eez, POI=POI,
                 fill_var=names(cmems_vars)[k], fill_lim=cmems_vars[[k]],
                 alpha_lim=NULL,#range(sqrt(abs(cmems_east$sla)), na.rm=T),
                 xlim=east_bbox[1:2], ylim=east_bbox[3:4],
                 title=paste0(i_time, ifelse(forecast, " FORECAST", "")),
                 out_dim=c(8, 4.5),
                 out_f=glue("{dir_figTemp}jc257_{ifelse(forecast, 'Fcst', 'Obs')}-",
                            "{ifelse(hull, 'hull', 'pk')}_{names(cmems_vars)[k]}_{i_time}.png"),
                 tracks=i_eddy,
                 darkLines=names(cmems_vars)[k]=="sla")
    }
    plot_cmems(df=i_df |> filter(!is.na(sla)), bathy=NULL, land=land, ccz=ccz, eez=eez, POI=POI,
               fill_var="uvDir", fill_lim=NULL,
               alpha_lim=c(0, max(cmems_vars$sla^2)),#range(sqrt(abs(cmems_east$sla)), na.rm=T),
               xlim=east_bbox[1:2], ylim=east_bbox[3:4],
               title=paste0(i_time, ifelse(forecast, " FORECAST", "")),
               out_dim=c(8, 4.5),
               out_f=glue("{dir_figTemp}jc257_{ifelse(forecast, 'Fcst', 'Obs')}-",
                          "{ifelse(hull, 'hull', 'pk')}_uvDir_{i_time}.png"),
               tracks=i_eddy,
               darkLines=T) 
  }
  gc()
}




img_path <- dir_figTemp
out_path <- "figs/cruise2024/"

fps <- 5
sets <- paste0("jc257_", ifelse(forecast, 'Fcst', 'Obs'), "-",
               ifelse(hull, 'hull', 'pk'), "_", c(names(cmems_vars), "uvDir"))

for(i in sets) {
  try({
    av_encode_video(dirf(img_path, i), 
                    glue("{out_path}{i}_{today()}.mp4"), 
                    framerate=fps)
  })
  gc()
}



if(FALSE) {
  f <- c(paste0("figs/cruise2024/forecast_2024-02-11_on_", today(), "_", c("EKE", "uv", "sla"), ".png"),
         paste0("figs/cruise2024/past_4wks.png"),
         paste0("figs/cruise2024/past_7days.png"),
         paste0("figs/cruise2024/SLA_trend_max.png"),
         paste0("figs/cruise2024/jc257_Obs-hull_", c("EKE", "uv", "sla"), "_", today(), ".mp4"),
         paste0("figs/cruise2024/jc257_Obs-pk_", c("EKE", "uv", "sla"), "_", today(), ".mp4"),
         paste0("D:/SMARTEX/figs/temp/jc257_Obs-hull_", c("EKE", "uv", "sla"), "_", today(), ".png"),
         paste0("D:/SMARTEX/figs/temp/jc257_Obs-pk_", c("EKE", "uv", "sla"), "_", today(), ".png"))
  walk(f, ~file.copy(.x, "P:/PHYSICS/SMARTEX/jc257/eddy_tracking/", overwrite=T))
}




# weekly change -----------------------------------------------------------

wk_seq <- seq(today()-3*7, today(), by=7)
p1 <- cmems_east |>
  filter(date %in% wk_seq) |>
  ggplot() + 
  geom_raster(aes(lon, lat, fill=sla)) +
  find_palette("sla", cmems_vars$sla) +
  geom_sf(data=eez, linewidth=0.75, colour="grey30", fill=NA) +
  geom_sf(data=ccz, linewidth=0.25, colour="grey30", fill=NA) +
  geom_sf(data=ccz |> filter(Contractor=="UKSRL"), linewidth=0.75, colour="grey30") +
  geom_point(data=POI, aes(lon, lat, shape=name), colour="grey30", size=3) +
  scale_shape_manual("", values=c(5,19,3,4), guide="none") +
  xlim(-120, -100) + ylim(east_bbox[[3]], east_bbox[[4]]) +
  facet_wrap(~date, ncol=1) +
  theme(axis.title=element_blank(),
        legend.position="bottom", 
        legend.key.height=unit(0.4, "cm"), 
        legend.key.width=unit(1.25, "cm"))
p2 <- cmems_east |>
  filter(date %in% wk_seq) |>
  ggplot() + 
  geom_raster(aes(lon, lat, fill=uv)) +
  find_palette("uv", cmems_vars$uv) +
  geom_sf(data=eez, linewidth=0.75, colour="grey30", fill=NA) +
  geom_sf(data=ccz, linewidth=0.25, colour="grey30", fill=NA) +
  geom_sf(data=ccz |> filter(Contractor=="UKSRL"), linewidth=0.75, colour="grey30") +
  geom_point(data=POI, aes(lon, lat, shape=name), colour="grey30", size=3) +
  scale_shape_manual("", values=c(5,19,3,4), guide="none") +
  xlim(-120, -100) + ylim(east_bbox[[3]], east_bbox[[4]]) +
  facet_wrap(~date, ncol=1) +
  theme(axis.title=element_blank(),
        legend.position="bottom", 
        legend.key.height=unit(0.4, "cm"), 
        legend.key.width=unit(1.25, "cm"))
ggpubr::ggarrange(p1, p2, ncol=2)
ggsave("figs/cruise2024/past_4wks.png", width=8, height=12, units="in")


wk_seq <- seq(today()-6, today(), by=2)
p1 <- cmems_east |>
  filter(date %in% wk_seq) |>
  ggplot() + 
  geom_raster(aes(lon, lat, fill=sla)) +
  find_palette("sla", cmems_vars$sla) +
  geom_sf(data=eez, linewidth=0.75, colour="grey30", fill=NA) +
  geom_sf(data=ccz, linewidth=0.25, colour="grey30", fill=NA) +
  geom_sf(data=ccz |> filter(Contractor=="UKSRL"), linewidth=0.75, colour="grey30") +
  geom_point(data=POI, aes(lon, lat, shape=name), colour="grey30", size=3) +
  scale_shape_manual("", values=c(5,19,3,4), guide="none") +
  xlim(-120, -100) + ylim(east_bbox[[3]], east_bbox[[4]]) +
  facet_wrap(~date, ncol=1) +
  theme(axis.title=element_blank(),
        legend.position="bottom", 
        legend.key.height=unit(0.4, "cm"), 
        legend.key.width=unit(1.25, "cm"))
p2 <- cmems_east |>
  filter(date %in% wk_seq) |>
  ggplot() + 
  geom_raster(aes(lon, lat, fill=uv)) +
  find_palette("uv", cmems_vars$uv) +
  geom_sf(data=eez, linewidth=0.75, colour="grey30", fill=NA) +
  geom_sf(data=ccz, linewidth=0.25, colour="grey30", fill=NA) +
  geom_sf(data=ccz |> filter(Contractor=="UKSRL"), linewidth=0.75, colour="grey30") +
  geom_point(data=POI, aes(lon, lat, shape=name), colour="grey30", size=3) +
  scale_shape_manual("", values=c(5,19,3,4), guide="none") +
  xlim(-120, -100) + ylim(east_bbox[[3]], east_bbox[[4]]) +
  facet_wrap(~date, ncol=1) +
  theme(axis.title=element_blank(),
        legend.position="bottom", 
        legend.key.height=unit(0.4, "cm"), 
        legend.key.width=unit(1.25, "cm"))
ggpubr::ggarrange(p1, p2, ncol=2)
ggsave("figs/cruise2024/past_7days.png", width=8, height=12, units="in")




# finding eddy centers ----------------------------------------------------

eddy_sf <- st_read("data/cruise2024/cruise2024_eddies.shp")
eddy_centroids <- eddy_sf |> 
  st_centroid() |>
  add_lonlat(drop_geom=T) |>
  rename(ctrd_lon=lon, ctrd_lat=lat) 
eddy_jc257 <- eddy_sf %>%
  mutate(rad=as.numeric(sqrt(st_area(.)/pi))) |>
  arrange(date) |>
  mutate(sla_max=NA_real_,
         peak_lon=NA_real_,
         peak_lat=NA_real_) |>
  left_join(eddy_centroids, by=c("id", "date")) |>
  as("SpatVector")
hull_ls <- hull90_ls <- vector("list", nrow(eddy_jc257))

for(i in 1:nrow(eddy_jc257)) {
  cmems_i <- glue("{dir_gpkg}interp/cmems_Obs_sla-interp_{eddy_jc257$date[i]}.tif") |>
    rast() |>
    crop(eddy_jc257[i,], mask=T)
  q95 <- global(cmems_i, fun=quantile, probs=0.95, na.rm=T)
  q90 <- global(cmems_i, fun=quantile, probs=0.9, na.rm=T)
  hull_ls[[i]] <- cmems_i |>
    as.points() |>
    st_as_sf() |>
    filter(sla.pred > q95[1,1]) |>
    summarise() |>
    st_convex_hull() |>
    mutate(date=eddy_jc257$date[i])
  hull90_ls[[i]] <- cmems_i |>
    as.points() |>
    st_as_sf() |>
    filter(sla.pred > q90[1,1]) |>
    summarise() |>
    st_convex_hull() |>
    mutate(date=eddy_jc257$date[i])
  pk_loc <- hull_ls[[i]] |> 
    st_centroid() |>
    st_coordinates()
  max_i <- where.max(cmems_i)
  eddy_jc257$sla_max[i] <- max_i[1,3]
  eddy_jc257$peak_lon[i] <- pk_loc[1,1]
  eddy_jc257$peak_lat[i] <- pk_loc[1,2]
}
eddy_jc257 <- eddy_jc257 |>
  as_tibble() 
eddy_jc257 |>
  st_as_sf(coords=c("peak_lon", "peak_lat"), crs=4326) |>
  write_sf("data/cruise2024/cruise2024_eddies_SSH_peak.gpkg")
eddy_jc257 |>
  st_as_sf(coords=c("ctrd_lon", "ctrd_lat"), crs=4326) %>%
  st_buffer(dist=.$rad) |>
  write_sf("data/cruise2024/cruise2024_eddies_radius.gpkg")
eddy_jc257 |>
  st_as_sf(coords=c("ctrd_lon", "ctrd_lat"), crs=4326) |>
  write_sf("data/cruise2024/cruise2024_eddies_centroids.gpkg")
hull_ls |> 
  reduce(bind_rows) |>
  write_sf(glue("{dir_gpkg}sla-q95-hull.gpkg"))
hull90_ls |> 
  reduce(bind_rows) |>
  write_sf(glue("{dir_gpkg}sla-q90-hull.gpkg"))










# quick autoreg -----------------------------------------------------------

eddy_jc257 <- st_read("data/cruise2024/cruise2024_eddies_SSH_peak.gpkg")
eddy_jc257_rad <- st_read("data/cruise2024/cruise2024_eddies_radius.gpkg")

# lm plots
ggplot(eddy_jc257, aes(date, sla_max)) + 
  geom_point() + 
  stat_smooth(data=eddy_jc257 |> filter(date >= "2023-12-31"), method="lm", formula=y~x+I(x^2)) +
  ylim(0.2, NA) + 
  scale_x_date(date_breaks="7 days") + 
  ylab("CMEMS SLA at eddy centroid (m)") + 
  theme_bw()
ggsave("figs/cruise2024/SLA_trend_max.png", width=8, height=4)
sla_lm <- lm(sla_max ~ yday + I(yday^2), 
             data=eddy_jc257 |> filter(date >= "2024-01-01") |> mutate(yday=yday(date)))
jc257_sla_pred <- tibble(date=seq(today(), ymd("2024-02-14"), by=1),
                         yday=yday(date)) %>%
  mutate(sla_pred=predict(sla_lm, newdata=.),
         sla_lo=predict(sla_lm, newdata=., interval="prediction")[,2],
         sla_hi=predict(sla_lm, newdata=., interval="prediction")[,3])
ggplot(eddy_jc257, aes(date, sla_max)) + 
  geom_point() + 
  stat_smooth(data=eddy_jc257 |> filter(date >= "2023-12-31"), method="lm", formula=y~x+I(x^2)) +
  geom_ribbon(data=jc257_sla_pred, aes(date, y=sla_pred, ymin=sla_lo, ymax=sla_hi), 
              alpha=0.25, colour=NA) +
  geom_point(data=jc257_sla_pred, aes(date, sla_pred), shape=1) +
  scale_x_date(date_breaks="7 days") + 
  ylim(0.2, NA) +
  ylab("CMEMS SLA at eddy centroid (m)") + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        axis.title.x=element_blank())
ggsave("figs/cruise2024/SLA_forecast_max.png", width=8, height=4)



# SSH Peak
eddy_jc257 <- st_read("data/cruise2024/cruise2024_eddies_SSH_peak.gpkg") 
hull_jc257 <- st_read(glue("{dir_gpkg}sla-q95-hull.gpkg")) %>%
  mutate(radius=as.numeric(sqrt(st_area(.)/pi)))
hull90_jc257 <- st_read(glue("{dir_gpkg}sla-q90-hull.gpkg")) %>%
  mutate(radius=as.numeric(sqrt(st_area(.)/pi)))
eddy_jc257 <- eddy_jc257 |>
  st_transform(8859) |>
  add_lonlat() |>
  mutate(dist=c(NA, st_distance(eddy_jc257[-1,], eddy_jc257[-nrow(eddy_jc257),], by_element=T)),
         dx=lon-lag(lon),
         dy=lat-lag(lat),
         theta=atan2(dy, dx)) 
summary(eddy_jc257)
pred_jc257 <- bind_rows(
  eddy_jc257 |> st_drop_geometry() |> slice_tail(n=1),
  tibble(id=1, 
         date=seq(today()+1, ymd("2024-02-11"), by=1),
         dist=median(eddy_jc257$dist, na.rm=T),
         dx=median(eddy_jc257$dx, na.rm=T),
         dy=median(eddy_jc257$dy, na.rm=T),
         theta=median(eddy_jc257$theta, na.rm=T),
         rad=median(eddy_jc257$rad)) 
)
for(i in 2:nrow(pred_jc257)) {
  pred_jc257$lon[i] <- pred_jc257$lon[i-1] + pred_jc257$dx[i]
  pred_jc257$lat[i] <- pred_jc257$lat[i-1] + pred_jc257$dy[i]
}
write_sf(pred_jc257 |> slice_tail(n=1) |> 
           st_as_sf(coords=c("lon", "lat"), crs=8859) |> st_transform(crs=4326), 
         glue("{dir_gpkg}forecasted_peak.gpkg"))



# Centroid + radius
eddy_jc257_ctrd <- st_read("data/cruise2024/cruise2024_eddies_centroids.gpkg")
eddy_jc257_ctrd <- eddy_jc257_ctrd |>
  arrange(date) |>
  st_transform(8859) |>
  add_lonlat() |>
  mutate(dist=c(NA, st_distance(eddy_jc257_ctrd[-1,], eddy_jc257_ctrd[-nrow(eddy_jc257_ctrd),], by_element=T)),
         dx=lon-lag(lon),
         dy=lat-lag(lat),
         theta=atan2(dy, dx)) 
summary(eddy_jc257_ctrd)
pred_jc257_ctrd <- bind_rows(
  eddy_jc257_ctrd |> st_drop_geometry() |> slice_tail(n=1),
  tibble(id=1, 
         date=seq(today()+1, ymd("2024-02-11"), by=1),
         dist=median(eddy_jc257_ctrd$dist, na.rm=T),
         dx=median(eddy_jc257_ctrd$dx, na.rm=T),
         dy=median(eddy_jc257_ctrd$dy, na.rm=T),
         theta=median(eddy_jc257_ctrd$theta, na.rm=T),
         rad=median(eddy_jc257_ctrd$rad)) 
)
for(i in 2:nrow(pred_jc257_ctrd)) {
  pred_jc257_ctrd$lon[i] <- pred_jc257_ctrd$lon[i-1] + pred_jc257_ctrd$dx[i]
  pred_jc257_ctrd$lat[i] <- pred_jc257_ctrd$lat[i-1] + pred_jc257_ctrd$dy[i]
}
pred_rad <- pred_jc257_ctrd |> 
  slice_tail(n=1) |>
  st_as_sf(coords=c("lon", "lat"), crs=8859) |>
  st_transform(4326) %>%
  st_buffer(dist=.$rad)
write_sf(pred_rad, glue("{dir_gpkg}forecasted_radius.gpkg"))


jc257_path <- st_read(glue("{dir_gpkg}jc257_proposed_path.shp")) |> filter(id < 3)

ggplot(eddy_jc257 |> st_transform(4326)) + 
  geom_raster(data=cmems_east |> filter(date==today()), aes(lon, lat, fill=sla)) +
  find_palette("sla", range(cmems_east$sla, na.rm=T)) +
  geom_sf(data=ccz_east) +
  geom_sf(data=eez, fill=NA, linewidth=1.5) +
  geom_sf() +
  geom_sf(data=pred_jc257 |> slice_tail(n=1) |> st_as_sf(coords=c("lon", "lat"), crs=8859) |>
            st_transform(4326) |> st_buffer(dist=median(hull_jc257$radius)), colour="red3", fill=NA) +
  geom_sf(data=pred_jc257 |> slice_tail(n=1) |> st_as_sf(coords=c("lon", "lat"), crs=8859) |>
            st_transform(4326) |> st_buffer(dist=median(hull90_jc257$radius)), colour="red3", fill=NA) +
  geom_sf(data=pred_rad, fill=NA, colour="red3") +
  geom_vline(xintercept=seq(-118, -105, by=1), colour="grey30", alpha=0.2) +
  geom_hline(yintercept=seq(9, 18, by=1), colour="grey30", alpha=0.2) +
  geom_sf(data=jc257_path, colour="#68228B", linewidth=1) +
  xlim(-118, -105) + ylim(9, 18) + 
  labs(title="2024-Feb-11 forecast based on median daily movement",
       subtitle=glue("Points: SSH peak; Inner circles: predicted top 5%, 10%; Altimetry: {today()}")) +
  theme(axis.title=element_blank())
ggsave(glue("figs/cruise2024/forecast_2024-02-11_on_{today()}_sla.png"), 
       width=8, height=5)
ggplot(eddy_jc257 |> st_transform(4326)) + 
  geom_raster(data=cmems_east |> filter(date==today()), aes(lon, lat, fill=uv)) +
  find_palette("uv", range(cmems_east$uv, na.rm=T)) +
  geom_sf(data=ccz_east) +
  geom_sf(data=eez, fill=NA, linewidth=1.5) +
  geom_sf() +
  geom_sf(data=pred_jc257 |> slice_tail(n=1) |> st_as_sf(coords=c("lon", "lat"), crs=8859) |>
            st_transform(4326) |> st_buffer(dist=median(hull_jc257$radius)), colour="red3", fill=NA) +
  geom_sf(data=pred_jc257 |> slice_tail(n=1) |> st_as_sf(coords=c("lon", "lat"), crs=8859) |>
            st_transform(4326) |> st_buffer(dist=median(hull90_jc257$radius)), colour="red3", fill=NA) +
  geom_sf(data=pred_rad, fill=NA, colour="red3") +
  geom_vline(xintercept=seq(-118, -105, by=1), colour="grey30", alpha=0.2) +
  geom_hline(yintercept=seq(9, 18, by=1), colour="grey30", alpha=0.2) +
  geom_sf(data=jc257_path, colour="#68228B", linewidth=1) +
  xlim(-118, -105) + ylim(9, 18) + 
  labs(title="2024-Feb-11 forecast based on median daily movement",
       subtitle=glue("Points: SSH peak; Inner circles: predicted top 5%, 10%; Altimetry: {today()}")) +
  theme(axis.title=element_blank())
ggsave(glue("figs/cruise2024/forecast_2024-02-11_on_{today()}_uv.png"), 
       width=8, height=5)
ggplot(eddy_jc257 |> st_transform(4326)) + 
  geom_raster(data=cmems_east |> filter(date==today()), aes(lon, lat, fill=EKE)) +
  find_palette("EKE", range(cmems_east$EKE, na.rm=T)) +
  geom_sf(data=ccz_east) +
  geom_sf(data=eez, fill=NA, linewidth=1.5) +
  geom_sf() +
  geom_sf(data=pred_jc257 |> slice_tail(n=1) |> st_as_sf(coords=c("lon", "lat"), crs=8859) |>
            st_transform(4326) |> st_buffer(dist=median(hull_jc257$radius)), colour="red3", fill=NA) +
  geom_sf(data=pred_jc257 |> slice_tail(n=1) |> st_as_sf(coords=c("lon", "lat"), crs=8859) |>
            st_transform(4326) |> st_buffer(dist=median(hull90_jc257$radius)), colour="red3", fill=NA) +
  geom_sf(data=pred_rad, fill=NA, colour="red3") +
  geom_vline(xintercept=seq(-118, -105, by=1), colour="grey30", alpha=0.2) +
  geom_hline(yintercept=seq(9, 18, by=1), colour="grey30", alpha=0.2) +
  geom_sf(data=jc257_path, colour="#68228B", linewidth=1) +
  xlim(-118, -105) + ylim(9, 18) + 
  labs(title="2024-Feb-11 forecast based on median daily movement",
       subtitle=glue("Points: SSH peak; Inner circles: predicted top 5%, 10%; Altimetry: {today()}")) +
  theme(axis.title=element_blank())
ggsave(glue("figs/cruise2024/forecast_2024-02-11_on_{today()}_EKE.png"), 
       width=8, height=5)






# depth profiles ----------------------------------------------------------

path_df <- bind_rows(
  bind_rows(
    read_csv("data/cruise2024/profile1_NOAA_bathymetry.csv") |>
      mutate(src="NOAA"),
    read_csv("data/cruise2024/profile1_GEBCO_bathymetry.csv") |>
      mutate(src="GEBCO")
  ) |> 
    mutate(LINE_ID="Northern") ,
  bind_rows(
    read_csv("data/cruise2024/profile2_NOAA_bathymetry.csv") |>
      mutate(src="NOAA"),
    read_csv("data/cruise2024/profile2_GEBCO_bathymetry.csv") |>
      mutate(src="GEBCO")
  ) |> 
    mutate(LINE_ID="Southern") 
) |>
  arrange(LINE_ID, X, src) |> 
  group_by(LINE_ID, X) |> 
  slice_tail(n=1) |> 
  ungroup()
 
read_csv("data/cruise2024/profile2_NOAA_bathymetry.csv") |>
  ggplot(aes(X, Z)) + 
  geom_path() + 
  # geom_rug(data=path_df |> filter(src=="NOAA"), colour="grey30", sides="b", linewidth=0.05) +
  scale_y_continuous(breaks=seq(-5000, -2500, by=500)) +
  labs(x="Longitude", y="Depth (m)",
       title="Depth profiles of proposed routes; NOAA multibeam (rug), else GEBCO") +
  theme_bw()
ggsave("figs/cruise2024/jc257_path_depth_profile_.png", width=12, height=5)




# 2023 eddy ---------------------------------------------------------------

cmems_2023 <- download_cmems(creds=readRDS("data/cmems_cred.rds"), 
                             bbox=east_bbox, 
                             dates=c(ymd("2022-12-01"), 
                                     ymd("2023-05-31")),
                             out_nc=glue("data/CMEMS/cmems_east_SMARTEX-eddy.nc"))

timesteps <- sort(unique(cmems_2023$date))
cmems_vars <- list("sla"=range(cmems_2023$sla, na.rm=T),
                   "uv"=range(cmems_2023$uv, na.rm=T))

for(i in seq_along(timesteps)) {
  i_time <- timesteps[i]
  i_df <- cmems_2023 |> 
    filter(date==i_time) |>
    mutate(id=row_number())
  
  for(k in seq_along(cmems_vars)) {
    plot_cmems(df=i_df |> filter(!is.na(sla)), bathy=NULL, land=land, ccz=ccz, eez=eez, POI=POI,
               fill_var=names(cmems_vars)[k], fill_lim=cmems_vars[[k]],
               alpha_lim=NULL,
               xlim=east_bbox[1:2], ylim=east_bbox[3:4],
               title=i_time,
               out_dim=c(10, 4.5),
               out_f=glue("{dir_figTemp}eddy2023_{names(cmems_vars)[k]}_{i_time}.png"),
               tracks=NULL,
               darkLines=names(cmems_vars)[k]=="sla")
  }
  gc()
}

img_path <- dir_figTemp
out_path <- "figs/cruise2024/"

fps <- 5
sets <- paste0("eddy2023_", c("sla", "uv"))

for(i in sets) {
  try({
    av_encode_video(dirf(img_path, i), 
                    glue("{out_path}{i}.mp4"), 
                    framerate=fps)
  })
  gc()
}





# mapped identifiers ------------------------------------------------------

eddy_dist <- read_csv("data/cruise2024/eddy_id_distances.csv") |>
  mutate(id=str_sub(ID_POINT, 1, 1),
         id_2=str_sub(ID_NEAR, 1, 1),
         date_1=ymd(str_sub(ID_POINT, 3, -1)),
         date_2=ymd(str_sub(ID_NEAR, 3, -1))) |>
  filter(id == id_2) |>
  arrange(id, date_1) |>
  group_by(id, date_1) |>
  slice_head(n=1) |>
  ungroup() |>
  select(id, date_1, date_2, DISTANCE)

eddy_dist <- expand_grid(id=1:4,
                    date=seq(ymd("2023-12-11"), ymd("2024-01-22"), by=7)) |>
  mutate(dist=c(54, 53, 25, 35, 49, 56, 50,
                28, 29, 42, 20, 33, 57, 19,
                62, 13, 111, 28, 74, 66, 62,
                23, 22, 50, 99, 17, 57, 54))
eddy_dist |>
  group_by(id) |>
  summarise(mn_wkly=mean(dist))

eddy_westward <- tibble(id=1:4,
                        total=c(306, 190, 331, 285)) |>
  mutate(mn_daily=total/as.numeric(ymd("2024-01-22") - ymd("2023-12-05"))) |>
  mutate(westward_byFeb09=mn_daily * as.numeric(ymd("2024-02-11") - ymd("2024-01-22")))














# hb model ----------------------------------------------------------------

d <- 1
ld_z <- readRDS(glue("data/cruise2024/aviso_ld_z_{d}d.rds"))
tracks_1day_jc257 <- eddy_sf %>%
  mutate(rad=as.numeric(sqrt(st_area(.)/pi))) |>
  st_centroid() |>
  mutate(track="jc257") |>
  add_lonlat() |>
  st_transform(8859) %>% # Equal Earth Asia-Pacific
  mutate(x=st_coordinates(.)[,1]/1e3,
         y=st_coordinates(.)[,2]/1e3,
         lrad=log(rad)) |>
  st_drop_geometry() |>
  arrange(track, date) |>
  group_by(track) |>
  mutate(x1=lag(x, d*1), y1=lag(y, d*1), lrad1=lag(lrad, d*1),
         x2=lag(x, d*2), y2=lag(y, d*2), lrad2=lag(lrad, d*2),
         x3=lag(x, d*3), y3=lag(y, d*3), lrad3=lag(lrad, d*3),
         x7=lag(x, d*7), y7=lag(y, d*7), lrad7=lag(lrad, d*7)) |>
  ungroup() |>
  rename(x0=x, y0=y) |>
  mutate(dx0=x0-x1, dy0=y0-y1,
         dx1=x1-x2, dy1=y1-y2,
         dx2=x2-x3, dy2=y2-y3,
         dx7=x0-x7, dy7=y0-y7,
         dist0=sqrt(dx0^2 + dy0^2),
         dist1=sqrt(dx1^2 + dy1^2),
         dist2=sqrt(dx2^2 + dy2^2),
         dist7=sqrt(dx7^2 + dy7^2),
         theta0=atan2(dy0, dx0),
         theta1=atan2(dy1, dx1),
         theta2=atan2(dy2, dx2),
         theta7=atan2(dy7, dx7)) %>%
  filter(complete.cases(.)) |>
  mutate(ldist0=log(dist0),
         ldist1=log(dist1),
         ldist2=log(dist2),
         ldist7=log(dist7)) |>
  mutate(altTheta0=if_else(theta0<0, theta0+pi, theta0-pi),
         altTheta1=if_else(theta1<0, theta1+pi, theta1-pi),
         altTheta2=if_else(theta2<0, theta2+pi, theta2-pi),
         altTheta7=if_else(theta7<0, theta7+pi, theta7-pi)) |>
  mutate(ldist0=(ldist0-ld_z$dist0[1])/ld_z$dist0[2],
         ldist1=(ldist1-ld_z$dist1[1])/ld_z$dist1[2],
         ldist2=(ldist2-ld_z$dist2[1])/ld_z$dist2[2],
         ldist7=(ldist7-ld_z$dist7[1])/ld_z$dist7[2]) 

library(brms)
brm_out <- readRDS(glue("data/cruise2024/aviso_{d}day-lag_brm.rds"))


get_xyr_post <- function(dat_df, preds, keep=c(""), recursive=F) {
  
  preds_ls <- list(p_ldist0=preds[,,1],
                   p_altTheta0=preds[,,2],
                   p_lrad=preds[,,3])
  preds_ls$p_theta0 <- preds_ls$p_altTheta0 + pi*-1*sign(preds_ls$p_altTheta0)
  preds_ls$p_dist0 <- pmin(exp(preds_ls$p_ldist0*ld_z$dist0[2] + ld_z$dist0[1]), 5e2)
  preds_ls$p_dx0 <- preds_ls$p_dist0 * cos(preds_ls$p_theta0)
  preds_ls$p_dy0 <- preds_ls$p_dist0 * sin(preds_ls$p_theta0)
  if(length(dim(preds_ls$p_dx0)) > 1) {
    preds_ls$p_x0 <- matrix(dat_df$x1, byrow=T,
                            nrow=nrow(preds_ls$p_dx0), ncol=ncol(preds_ls$p_dx0)) + preds_ls$p_dx0
    preds_ls$p_y0 <- matrix(dat_df$y1, byrow=T,
                            nrow=nrow(preds_ls$p_dy0), ncol=ncol(preds_ls$p_dy0)) + preds_ls$p_dy0
  } else {
    preds_ls$p_x0 <- dat_df$x1 + preds_ls$p_dx0
    preds_ls$p_y0 <- dat_df$y1 + preds_ls$p_dy0
  }
  
  preds_ls$p_rad <- exp(preds_ls$p_lrad)
  
  if(recursive) {
    pred_df <- expand_grid(date=unique(dat_df$date),
                           iter1=1:dim(preds)[2],
                           iter=1:dim(preds)[1]) |>
      mutate(p_x0=c(preds_ls$p_x0),
             p_y0=c(preds_ls$p_y0),
             p_rad=c(preds_ls$p_rad),
             p_ldist0=c(preds_ls$p_ldist0),
             p_lrad=c(preds_ls$p_lrad),
             p_altTheta0=c(preds_ls$p_altTheta0))|>
      left_join(dat_df |> 
                  select(any_of(c(keep))) |> 
                  mutate(iter1=row_number()),
                by=c("date", "iter1")) |>
      mutate(iter=row_number()) |> select(-iter1)
  } else {
    pred_df <- expand_grid(date=unique(dat_df$date),
                           iter=1:dim(preds)[1]) |>
      mutate(p_x0=c(preds_ls$p_x0),
             p_y0=c(preds_ls$p_y0),
             p_rad=c(preds_ls$p_rad),
             p_ldist0=c(preds_ls$p_ldist0),
             p_lrad=c(preds_ls$p_lrad),
             p_altTheta0=c(preds_ls$p_altTheta0))|>
      left_join(dat_df |> select(any_of(keep)))
  }
  
  return(pred_df)
}


pred_t5 <- tracks_1day_jc257 %>%
  get_xyr_post(., posterior_epred(brm_out, newdata=.), 
               keep=c("date", "track", "x0", "y0", "rad"))

iter_samp <- sample(1:max(pred_t5$iter), 100)

pred_t5_sf <- pred_t5 |>
  filter(iter %in% iter_samp) |>
  mutate(p_x0=p_x0*1e3, p_y0=p_y0*1e3) |>
  st_as_sf(coords=c("p_x0", "p_y0"), crs=8859) |>
  st_transform(4326) %>%
  st_buffer(dist=.$p_rad)

ggplot(pred_t5_sf) + 
  geom_sf(fill="grey30", colour=NA, alpha=0.01) + 
  geom_sf(data=eddy_sf |> st_centroid()) + 
  geom_sf(data=ccz, fill=NA, colour="black")

jc257_timesteps <- seq(max(tracks_1day_jc257$date)+1, ymd("2024-02-14"), by=1)

pred_ls <- vector("list", length(jc257_timesteps))
for(i in seq_along(jc257_timesteps)) {
  i_date <- jc257_timesteps[i]
  if(i == 1) {
    last_df <- tracks_1day_jc257 |>
      filter(date==i_date-1) |> mutate(iter=1)
  } else {
    iters <- sample(unique(pred_ls[[i-1]]$iter), 500)
    last_df <- pred_ls[[i-1]] |> filter(iter %in% iters) |>
      rename(x0=p_x0, y0=p_y0, lrad=p_lrad, ldist0=p_ldist0, altTheta0=p_altTheta0)
  }

  pred_ls[[i]] <- last_df |>
    mutate(date=i_date) |>
    select(track, iter, date,
           x0, y0, lrad, ldist0, altTheta0,
           x1, y1, lrad1, ldist1, altTheta1) |>
    rename(x2=x1, y2=y1, ldist2=ldist1, lrad2=lrad1, altTheta2=altTheta1,
           x1=x0, y1=y0, ldist1=ldist0, lrad1=lrad, altTheta1=altTheta0) %>%
    get_xyr_post(posterior_epred(brm_out, newdata=.), 
                 keep=c("date", "track",
                        "x1", "y1", "ldist1", "lrad1", "altTheta1",
                        "x2", "y2", "ldist2", "lrad2", "altTheta2"),
                 recursive=T)
  gc()
}

n_iters <- 3e3
forecast_df <- bind_rows(
  pred_t5 |> sample_n(n_iters), 
  map_dfr(pred_ls, ~.x |> sample_n(n_iters)) |>
    mutate(p_theta0=p_altTheta0 + pi*if_else(p_altTheta0>0, -1, 1),
           p_dist0=pmin(exp(p_ldist0*ld_z$dist0[2] + ld_z$dist0[1]), 5e2),
           p_dx0=p_dist0 * cos(p_theta0),
           p_dy0=p_dist0 * sin(p_theta0)) 
  )
ggplot(forecast_df, aes(p_x0, p_y0)) + geom_point(alpha=0.1, shape=1)

pred_sf <- forecast_df |>
  filter(date > today()) |>
  group_by(date) |>
  sample_n(pmin(n(), 100)) |>
  ungroup() |>
  mutate(p_x0=p_x0*1e3, p_y0=p_y0*1e3) |>
  st_as_sf(coords=c("p_x0", "p_y0"), crs=8859) |>
  st_transform(4326)  %>%
  st_buffer(dist=.$p_rad)
obs_sf <- pred_t5 |> filter(iter==1) |>
  mutate(x0=x0*1e3, y0=y0*1e3) |>
  st_as_sf(coords=c("x0", "y0"), crs=8859) |>
  st_transform(4326)  %>%
  st_buffer(dist=.$rad)
pred_sf |>
  ggplot() +
  geom_sf(fill="grey30", alpha=0.04, colour=NA) +
  geom_sf(data=pred_sf |> st_centroid(), alpha=0.04) +
  geom_sf(data=ccz, colour="black") +
  geom_sf(data=obs_sf, fill=NA, colour="black") + 
  scale_colour_brewer(type="qual", palette=2)




pred_t5 |>
  ggplot(aes(p_x0, x0)) + geom_point() + geom_abline()
pred_t5 |>
  ggplot(aes(p_y0, y0)) + geom_point() + geom_abline()

pred_t5 |>
  ggplot(aes(p_ldist0, ldist0)) + geom_point() + geom_abline()
pred_t5 |>
  ggplot(aes(p_dist0, dist0)) + geom_point() + geom_abline()

pred_t5 |>
  ggplot(aes(p_lrad, lrad)) + geom_point() + geom_abline()
pred_t5 |>
  ggplot(aes(p_rad, rad)) + geom_point() + geom_abline()
pred_t5 |>
  ggplot(aes(p_altTheta0, altTheta0)) + geom_point() + geom_abline()

pred_t5 |>
  ggplot(aes((p_ldist0 - ldist0)/abs(ldist0))) + geom_histogram()




eddy_jc257 %>%
  st_buffer(dist=.$rad) |>
  st_transform(4326) |>
  ggplot() +
  geom_sf(data=bathy, aes(fill=elev), colour=NA) +
  scale_fill_continuous_divergingx(h1=0, h2=NA, h3=0, 
                                   c1=0, c2=NA, c3=0, 
                                   l1=0, l2=90, l3=0, 
                                   p1=1.3, p2=NA, p3=1.3, p4=0.75, 
                                   cmax1=NA, cmax2=NA, guide="none") +
  geom_sf(data=ccz, colour="cadetblue") +
  geom_sf(fill=NA, alpha=0.5, aes(colour=factor(id))) + 
  scale_colour_brewer(type="qual", palette=2)

eddy_jc257 %>%
  st_buffer(dist=.$rad) |>
  ggplot() +
  geom_sf(data=bathy, aes(fill=elev), colour=NA) +
  scale_fill_continuous_divergingx(h1=0, h2=NA, h3=0, 
                                   c1=0, c2=NA, c3=0, 
                                   l1=0, l2=90, l3=0, 
                                   p1=1.3, p2=NA, p3=1.3, p4=0.75, 
                                   cmax1=NA, cmax2=NA, guide="none") +
  geom_sf(data=ccz, colour="cadetblue") +
  geom_sf(fill=NA, colour="red", alpha=0.5) + 
  facet_wrap(~id)




