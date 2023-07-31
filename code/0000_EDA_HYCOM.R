
rm(list=ls()); gc()

library(tidyverse)
library(lubridate)
library(ncdf4)
library(oce)
library(sf)
library(gganimate)
theme_set(theme_classic())


ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
ccz_bbox <- st_bbox(ccz) |> as.list()
ccz_crossSections <- tibble(y0=ccz_bbox$ymin,
                           y1=ccz_bbox$ymax,
                           x0=seq(ccz_bbox$xmin, ccz_bbox$xmax, length.out=5)) |>
  mutate(x1=x0) |>
  rowwise() |>
  mutate(geom=c(x0, x1, y0, y1) |>
           matrix(nrow=2) |>
           st_linestring() |>
           st_sfc(crs=4326) |>
           c()) |>
  ungroup() |>
  st_as_sf()
  
dates <- c("2022-06-01", "2022-06-03")
dates_s <- difftime(as_datetime(dates), 
                    as_datetime("2000-01-01 00:00:00"), 
                    units="hours") |>
  as.numeric()
thredds_url <- "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z"
hy_nc <- nc_open(thredds_url)
sub_i <- map(1:nrow(ccz_crossSections),
             ~list(
               lon=abs(ncvar_get(hy_nc, "lon") - (360 + ccz_crossSections$x0[.x])) |>
                 which.min(),
               lat=between(ncvar_get(hy_nc, "lat"), ccz_crossSections$y0[.x], ccz_crossSections$y1[.x]) |> 
                 which() |> range(),
               time=between(ncvar_get(hy_nc, "time"), dates_s[1], dates_s[2]) |> which() |> range(),
               depth=c(0, 39)
             ) |>
               map(~paste0("[", min(.x), ":", max(.x), "]"))
)
nc_close(hy_nc)

suffix <- map(sub_i,
              ~paste0("?",
                      "depth", .x$depth, ",",
                      "lon", .x$lon, ",",
                      "lat", .x$lat, ",",
                      "time", .x$time, ",",
                      "water_u", .x$time, .x$depth, .x$lat, .x$lon, ",",
                      "water_v", .x$time, .x$depth, .x$lat, .x$lon, ",",
                      "water_u_bottom", .x$time, .x$lat, .x$lon, ",",
                      "water_v_bottom", .x$time, .x$lat, .x$lon))
hy_ls <- hy_curl <- hy_df <- map(seq_along(suffix), ~list())
for(i in seq_along(suffix)) {
  nc_url <- paste0(thredds_url, suffix[[i]])
  hy_nc <- nc_open(nc_url)
  hy_ls[[i]] <- list(
    depth=c(ncvar_get(hy_nc, "depth")),
    time=c(ncvar_get(hy_nc, "time")),
    lon=c(ncvar_get(hy_nc, "lon")),
    lat=c(ncvar_get(hy_nc, "lat")),
    u=ncvar_get(hy_nc, "water_u"),
    v=ncvar_get(hy_nc, "water_v")
  )
  nc_close(hy_nc)
  
  hy_df[[i]] <- expand_grid(time=hy_ls[[i]]$time,
                            depth=hy_ls[[i]]$depth,
                            lat=hy_ls[[i]]$lat,
                            lon=hy_ls[[i]]$lon) |>
    arrange(lon, lat, depth, time) |>
    mutate(id=as.numeric(as.factor(paste(lon, lat))),
           time=as_datetime("2000-01-01 00:00:00") + time*60*60,
           u=c(hy_ls[[i]]$u),
           v=c(hy_ls[[i]]$v),
           curl=) |>
    mutate(uv=sqrt(u^2 + v^2))
  
}


hy_curl <- map_dfr(seq_along(hy_ls$time),
                   ~expand_grid(lat=hy_ls$lat, 
                                lon=hy_ls$lon) |>
                     mutate(curl_sur=c(curl(hy_ls$u[,,.x], 
                                            hy_ls$v[,,.x], 
                                            hy_ls$lon, 
                                            hy_ls$lat, 
                                            geographical=T)$curl),
                            curl_btm=c(curl(hy_ls$u_bottom[,,.x], 
                                            hy_ls$v_bottom[,,.x], 
                                            hy_ls$lon, 
                                            hy_ls$lat, 
                                            geographical=T)$curl),
                            time=hy_ls$time[.x])) |>
  arrange(lon, lat, time) |>
  mutate(id=as.numeric(as.factor(paste(lon, lat))),
         time=as_datetime("2000-01-01 00:00:00") + time*60*60) |>
  mutate(across(starts_with("curl"), ~.x/sd(.x, na.rm=T), .names="{.col}_std"))

hy_df <- expand_grid(time=hy_ls$time,
                      # depth=hy_ls$depth,
                      lat=hy_ls$lat,
                      lon=hy_ls$lon) |>
  arrange(lon, lat, time) |>
  mutate(id=as.numeric(as.factor(paste(lon, lat))),
         time=as_datetime("2000-01-01 00:00:00") + time*60*60,
         u_sur=c(hy_ls$u),
         v_sur=c(hy_ls$v),
         u_btm=c(hy_ls$u_btm),
         v_btm=c(hy_ls$v_btm)) |>
  mutate(uv_sur=sqrt(u_sur^2 + v_sur^2),
         uv_btm=sqrt(u_btm^2 + v_btm^2))

anim <- hy_curl |>
  ggplot(aes(lon, lat, fill=curl_sur_std, group=id)) + 
  geom_tile() + 
  scale_fill_distiller(type="div", palette="RdBu") +
  transition_states(time, transition_length=0, state_length=1) + 
  labs(title="{closest_state}")
anim_save("figs/anim/00_HYCOM_test_curl_sur.gif", anim, 
          nframes=n_distinct(hy_curl$time), fps=10, 
          width=6, height=4, units="in", res=300)

anim <- hy_curl |>
  ggplot(aes(lon, lat, fill=curl_btm_std, group=id)) + 
  geom_tile() + 
  scale_fill_distiller(type="div", palette="RdBu") +
  transition_states(time, transition_length=0, state_length=1) + 
  labs(title="{closest_state}")
anim_save("figs/anim/00_HYCOM_test_curl_btm.gif", anim, 
          nframes=n_distinct(hy_curl$time), fps=10, 
          width=6, height=4, units="in", res=300)

anim <- hy_df |>
  ggplot(aes(lon, lat, colour=uv_sur, group=id)) + 
  geom_point() + 
  scale_colour_viridis_c() +
  transition_states(time, transition_length=0, state_length=1) + 
  labs(title="{closest_state}")
anim_save("figs/anim/00_HYCOM_test_uv_sur.gif", anim, 
          nframes=n_distinct(hy_df$time), fps=10, 
          width=6, height=4, units="in", res=300)

anim <- hy_df |>
  ggplot(aes(lon, lat, colour=uv_btm, group=id)) + 
  geom_point() + 
  scale_colour_viridis_c() +
  transition_states(time, transition_length=0, state_length=1) + 
  labs(title="{closest_state}")
anim_save("figs/anim/00_HYCOM_test_uv_btm.gif", anim, 
          nframes=n_distinct(hy_df$time), fps=10, 
          width=6, height=4, units="in", res=300)

