



# setup -------------------------------------------------------------------

rm(list=ls()); gc()

library(tidyverse)
library(glue)
library(ncdf4)
library(oce)
library(sf)
library(magick)
library(colorspace)
library(sevcheck)
theme_set(theme_bw() + 
            theme(panel.grid=element_blank(),
                  legend.position="bottom",
                  legend.key.width=unit(2, "cm"),
                  legend.key.height=unit(0.25, "cm")))
cmr <- readRDS("data/cmr_cmaps.RDS")
walk(dir("code/fn", full.names=T), source)

E_dir <- "E:/Projects/SMARTEX/SMARTEX/"

ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
ccz_bbox <- st_bbox(ccz) |> as.list()
ccz_bbox$xmin <- -130
ccz_bbox$xmax <- -90

ccz_lat_slice <- tibble(y0=ccz_bbox$ymin,
                        y1=ccz_bbox$ymax,
                        x0=seq(ccz_bbox$xmin, ccz_bbox$xmax, length.out=5)) |>
  mutate(x1=x0) |>
  rowwise() |>
  mutate(geom=c(x0, x1, y0, y1) |> matrix(nrow=2) |>
           st_linestring() |> st_sfc(crs=4326) |> c()) |>
  ungroup() |> st_as_sf()
ccz_lon_slice <- tibble(x0=ccz_bbox$xmin,
                        x1=ccz_bbox$xmax,
                        y0=seq(ccz_bbox$ymin, ccz_bbox$ymax, length.out=5)) |>
  mutate(y1=y0) |>
  rowwise() |>
  mutate(geom=c(x0, x1, y0, y1) |> matrix(nrow=2) |>
           st_linestring() |> st_sfc(crs=4326) |> c()) |>
  ungroup() |> st_as_sf()
ccz_dep_slice <- tibble(x0=ccz_bbox$xmin, 
                        x1=ccz_bbox$xmax,
                        y0=ccz_bbox$ymin,
                        y1=ccz_bbox$ymax,
                        z=c(0,100,500,1000,1500,2000,3000,4000))


date_df <- tibble(start=ymd("2022-10-01") + seq(0, 303, by=1),
                  end=ymd("2022-10-01") + seq(0, 303, by=1))



# hycom indexes -----------------------------------------------------------

thredds_url <- "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z"
hy_nc <- nc_open(thredds_url)
hy_i <- list(lon=ncvar_get(hy_nc, "lon"),
             lat=ncvar_get(hy_nc, "lat"),
             time=ncvar_get(hy_nc, "time"),
             depth=ncvar_get(hy_nc, "depth"))
nc_close(hy_nc)



# get data ----------------------------------------------------------------

for(j in 194:nrow(date_df)) {
  d <- c(date_df[j,])
  dates <- c(paste0(year(d[[1]]), "-", month(d[[1]]), "-", day(d[[1]]), " 00:00:00"),
             paste0(year(d[[2]]), "-", month(d[[2]]), "-", day(d[[2]]), " 23:59:59"))
  dates_h <- as.numeric(difftime(as_datetime(dates),
                                 as_datetime("2000-01-01 00:00:00"),
                                 units="hours"))
  cat("Starting", as.character(d$start), "at", as.character(Sys.time()), "\n")

# * lat slices ------------------------------------------------------------

  # lat_i <- map(1:nrow(ccz_lat_slice),
  #              ~list(
  #                lon=abs(hy_i$lon - (360 + ccz_lat_slice$x0[.x])) |> which.min(),
  #                lat=between(hy_i$lat, ccz_lat_slice$y0[.x], ccz_lat_slice$y1[.x]) |>
  #                  which() |> range(),
  #                time=between(hy_i$time, dates_h[1], dates_h[2]) |> which() |> range(),
  #                depth=c(1, 40)) |>
  #                map(~paste0("[", min(.x)-1, ":", max(.x)-1, "]")))
  # lat_suf <- map(lat_i,
  #                ~paste0("?",
  #                        "depth", .x$depth, ",",
  #                        "lon", .x$lon, ",",
  #                        "lat", .x$lat, ",",
  #                        "time", .x$time, ",",
  #                        "water_u", .x$time, .x$depth, .x$lat, .x$lon, ",",
  #                        "water_v", .x$time, .x$depth, .x$lat, .x$lon))
  # 
  # lat_ls <- lat_df.ls <- map(seq_along(lat_suf), ~list())
  # for(i in seq_along(lat_suf)) {
  #   nc_url <- paste0(thredds_url, lat_suf[[i]])
  #   hy_nc <- nc_open(nc_url)
  #   lat_ls[[i]] <- list(
  #     depth=c(ncvar_get(hy_nc, "depth")),
  #     time=c(ncvar_get(hy_nc, "time")),
  #     lon=c(ncvar_get(hy_nc, "lon")),
  #     lat=c(ncvar_get(hy_nc, "lat")),
  #     u=ncvar_get(hy_nc, "water_u"),
  #     v=ncvar_get(hy_nc, "water_v")
  #   )
  #   nc_close(hy_nc)
  # 
  #   depth_layers <- tibble(depth=lat_ls[[i]]$depth) |>
  #     mutate(depth0=depth - (depth-lag(depth))/2,
  #            depth1=depth + (lead(depth)-depth)/2) |>
  #     mutate(depth0=if_else(is.na(depth0), min(depth), depth0),
  #            depth1=if_else(is.na(depth1), max(depth), depth1))
  #   lat_diff <- mean(lat_ls[[i]]$lat - lag(lat_ls[[i]]$lat), na.rm=T)
  #   lat_df.ls[[i]] <- expand_grid(time=lat_ls[[i]]$time,
  #                                depth=lat_ls[[i]]$depth,
  #                                lat=lat_ls[[i]]$lat,
  #                                lon=lat_ls[[i]]$lon) |>
  #     arrange(time, depth, lat, lon) |>
  #     mutate(id=as.numeric(as.factor(paste(lon, lat))),
  #            time=as_datetime("2000-01-01 00:00:00") + time*60*60,
  #            u=c(lat_ls[[i]]$u),
  #            v=c(lat_ls[[i]]$v)) |>
  #     mutate(uv=sqrt(u^2 + v^2)) |>
  #     mutate(lat_diff=lat_diff) |>
  #     left_join(depth_layers)
  #   gc()
  # }
  # latSlice_df <- do.call('bind_rows', lat_df.ls) |>
  #   mutate(lon=round(lon, 2),
  #          uvDir=atan2(v,u))
  # saveRDS(latSlice_df, glue("data/HYCOM/lat_slices_{d[[1]]}.rds"))
  # rm(lat_df.ls); rm(lat_ls); rm(latSlice_df)
  # gc()
  
  

# * lon slices ------------------------------------------------------------
  
  # lon_i <- map(1:nrow(ccz_lon_slice),
  #              ~list(
  #                lon=between(hy_i$lon, 360+ccz_lon_slice$x0[.x], 360+ccz_lon_slice$x1[.x]) |>
  #                  which() |> range(),
  #                lat=abs(hy_i$lat - ccz_lon_slice$y0[.x]) |> which.min(),
  #                time=between(hy_i$time, dates_h[1], dates_h[2]) |> which() |> range(),
  #                depth=c(1, 40)) |>
  #                map(~paste0("[", min(.x)-1, ":", max(.x)-1, "]")))
  # lon_suf <- map(lon_i,
  #                ~paste0("?",
  #                        "depth", .x$depth, ",",
  #                        "lon", .x$lon, ",",
  #                        "lat", .x$lat, ",",
  #                        "time", .x$time, ",",
  #                        "water_u", .x$time, .x$depth, .x$lat, .x$lon, ",",
  #                        "water_v", .x$time, .x$depth, .x$lat, .x$lon))
  # 
  # lon_ls <- lon_df.ls <- map(seq_along(lon_suf), ~list())
  # for(i in seq_along(lon_suf)) {
  #   nc_url <- paste0(thredds_url, lon_suf[[i]])
  #   hy_nc <- nc_open(nc_url)
  #   lon_ls[[i]] <- list(
  #     depth=c(ncvar_get(hy_nc, "depth")),
  #     time=c(ncvar_get(hy_nc, "time")),
  #     lon=c(ncvar_get(hy_nc, "lon")),
  #     lat=c(ncvar_get(hy_nc, "lat")),
  #     u=ncvar_get(hy_nc, "water_u"),
  #     v=ncvar_get(hy_nc, "water_v")
  #   )
  #   nc_close(hy_nc)
  # 
  #   depth_layers <- tibble(depth=lon_ls[[i]]$depth) |>
  #     mutate(depth0=depth - (depth-lag(depth))/2,
  #            depth1=depth + (lead(depth)-depth)/2) |>
  #     mutate(depth0=if_else(is.na(depth0), min(depth), depth0),
  #            depth1=if_else(is.na(depth1), max(depth), depth1))
  #   lon_diff <- mean(lon_ls[[i]]$lon - lag(lon_ls[[i]]$lon), na.rm=T)
  #   lon_df.ls[[i]] <- expand_grid(time=lon_ls[[i]]$time,
  #                                depth=lon_ls[[i]]$depth,
  #                                lat=lon_ls[[i]]$lat,
  #                                lon=lon_ls[[i]]$lon) |>
  #     arrange(time, depth, lat, lon) |>
  #     mutate(id=as.numeric(as.factor(paste(lon, lat))),
  #            time=as_datetime("2000-01-01 00:00:00") + time*60*60,
  #            u=c(lon_ls[[i]]$u),
  #            v=c(lon_ls[[i]]$v)) |>
  #     mutate(uv=sqrt(u^2 + v^2)) |>
  #     mutate(lon_diff=lon_diff) |>
  #     left_join(depth_layers)
  #   gc()
  # }
  # lonSlice_df <- do.call('bind_rows', lon_df.ls) |>
  #   mutate(lon=round(lon, 2),
  #          uvDir=atan2(v,u)) |>
  #   mutate(latFac=factor(lat, levels=rev(unique(lat))))
  # saveRDS(lonSlice_df, glue("data/HYCOM/lon_slices_{d[[1]]}.rds"))
  # rm(lon_df.ls); rm(lon_ls); rm(lonSlice_df)
  # gc()
  
  
# * depth slices ----------------------------------------------------------
  
  dep_i <- map(1:nrow(ccz_dep_slice),
               ~list(
                 lon=between(hy_i$lon, 360+ccz_dep_slice$x0[.x], 360+ccz_dep_slice$x1[.x]) |>
                   which() |> range(),
                 lat=between(hy_i$lat, ccz_dep_slice$y0[.x], ccz_dep_slice$y1[.x]) |>
                   which() |> range(),
                 time=between(hy_i$time, dates_h[1], dates_h[2]) |> which() |> range(),
                 depth=abs(hy_i$depth - ccz_dep_slice$z[.x]) |> which.min()) |>
                 map(~paste0("[", min(.x)-1, ":", max(.x)-1, "]")))
  dep_suf <- map(dep_i,
                 ~paste0("?",
                         "depth", .x$depth, ",",
                         "lon", .x$lon, ",",
                         "lat", .x$lat, ",",
                         "time", .x$time, ",",
                         "water_u", .x$time, .x$depth, .x$lat, .x$lon, ",",
                         "water_v", .x$time, .x$depth, .x$lat, .x$lon))

  dep_ls <- dep_df.ls <- map(seq_along(dep_suf), ~list())
  for(i in seq_along(dep_suf)) {
    nc_url <- paste0(thredds_url, dep_suf[[i]])
    hy_nc <- nc_open(nc_url)
    dep_ls[[i]] <- list(
      depth=c(ncvar_get(hy_nc, "depth")),
      time=c(ncvar_get(hy_nc, "time")),
      lon=c(ncvar_get(hy_nc, "lon")),
      lat=c(ncvar_get(hy_nc, "lat")),
      u=ncvar_get(hy_nc, "water_u"),
      v=ncvar_get(hy_nc, "water_v")
    )
    nc_close(hy_nc)

    dep_df.ls[[i]] <- expand_grid(time=dep_ls[[i]]$time,
                                  depth=dep_ls[[i]]$depth,
                                  lat=dep_ls[[i]]$lat,
                                  lon=dep_ls[[i]]$lon) |>
      arrange(time, depth, lat, lon) |>
      mutate(id=as.numeric(as.factor(paste(lon, lat))),
             time=as_datetime("2000-01-01 00:00:00") + time*60*60,
             u=c(dep_ls[[i]]$u),
             v=c(dep_ls[[i]]$v)) |>
      mutate(uv=sqrt(u^2 + v^2))
    gc()
  }
  depSlice_df <- do.call('bind_rows', dep_df.ls) |>
    mutate(uvDir=atan2(v,u))
  saveRDS(depSlice_df, glue("{E_dir}data/HYCOM/dep_slices_{d[[1]]}.rds"))
  rm(dep_df.ls); rm(dep_ls); rm(depSlice_df)
  gc()
}





# animations --------------------------------------------------------------

# * lat slices ------------------------------------------------------------

# latSlice_df <- map_dfr(dir("data/HYCOM", "lat_slices_", full.names=T), readRDS) |>
#   group_by(depth) |>
#   mutate(u_std=c(scale(u, center=F)),
#          v_std=c(scale(v, center=F)),
#          uv_std=c(scale(uv, center=F))) |>
#   ungroup()
# gc()
# 
# lims <- latSlice_df |>
#   select(u, v, uv, u_std, v_std, uv_std) |>
#   reframe(across(everything(), ~range(.x, na.rm=T)))
# 
# timesteps <- unique(latSlice_df$time)
# 
# for(i in seq_along(timesteps)) {
#   df_i <- latSlice_df |> 
#     filter(time==timesteps[i]) |>
#     mutate(lat0=lat-lat_diff,
#            lat1=lat+lat_diff)
#   
#   # * * u -------------------------------------------------------------------
#   p <- df_i |>
#     ggplot(aes(ymin=lat0, ymax=lat1, xmin=depth0, xmax=depth1, fill=u)) + 
#     geom_rect() + 
#     scale_fill_continuous_divergingx(palette="RdBu", na.value="grey85",
#                                      mid=0, limits=lims$u) + 
#     facet_grid(.~lon) + 
#     labs(title=timesteps[i], x="Depth (m)", y="Latitude") + 
#     scale_x_continuous(breaks=c(0, 2000, 4000))
#   ggsave(glue("figs/anim/temp/HYCOM_lat_slices_u-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=6, units="in", dpi=300)
# 
#   # * * v -------------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(ymin=lat0, ymax=lat1, xmin=depth0, xmax=depth1, fill=v)) + 
#     geom_rect() + 
#     scale_fill_continuous_divergingx(palette="PRGn", na.value="grey85", 
#                                      mid=0, limits=lims$v) + 
#     facet_grid(.~lon) + 
#     labs(title=timesteps[i], x="Depth (m)", y="Latitude") + 
#     scale_x_continuous(breaks=c(0, 2000, 4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lat_slices_v-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=6, units="in", dpi=300)
# 
#   # * * uv ------------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(ymin=lat0, ymax=lat1, xmin=depth0, xmax=depth1, fill=uv)) + 
#     geom_rect() + 
#     scale_fill_viridis_c(na.value="grey85", limits=lims$uv) + 
#     facet_grid(.~lon) + 
#     labs(title=timesteps[i], x="Depth (m)", y="Latitude") + 
#     scale_x_continuous(breaks=c(0, 2000, 4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lat_slices_uv-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=6, units="in", dpi=300)
# 
#   # * * uvDir ---------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(ymin=lat0, ymax=lat1, xmin=depth0, xmax=depth1, fill=uvDir)) + 
#     geom_rect() + 
#     scale_fill_gradientn(colours=cmr$seasons, limits=c(-pi, pi),
#                          breaks=c(-pi, -pi/2, 0, pi/2, pi), 
#                          labels=c("W", "S", "E", "N", "W")) + 
#     facet_grid(.~lon) + 
#     labs(title=timesteps[i], x="Depth (m)", y="Latitude") + 
#     scale_x_continuous(breaks=c(0, 2000, 4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lat_slices_uvDir-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=6, units="in", dpi=300)
#   
#   # * * u_std ---------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(ymin=lat0, ymax=lat1, xmin=depth0, xmax=depth1, fill=u_std)) + 
#     geom_rect() + 
#     scale_fill_continuous_divergingx(palette="RdBu", na.value="grey85", 
#                                      mid=0, limits=lims$u_std) + 
#     facet_grid(.~lon) + 
#     labs(title=timesteps[i], x="Depth (m)", y="Latitude") + 
#     scale_x_continuous(breaks=c(0, 2000, 4000))
#   ggsave(glue("figs/anim/temp/HYCOM_lat_slices_u_std-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=6, units="in", dpi=300)
#   
# 
#   # * * v_std ---------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(ymin=lat0, ymax=lat1, xmin=depth0, xmax=depth1, fill=v_std)) + 
#     geom_rect() + 
#     scale_fill_continuous_divergingx(palette="PRGn", na.value="grey85",
#                                      mid=0, limits=lims$v_std) + 
#     facet_grid(.~lon) + 
#     labs(title=timesteps[i], x="Depth (m)", y="Latitude") + 
#     scale_x_continuous(breaks=c(0, 2000, 4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lat_slices_v_std-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=6, units="in", dpi=300)
# 
#   # * * uv_std --------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(ymin=lat0, ymax=lat1, xmin=depth0, xmax=depth1, fill=uv_std)) + 
#     geom_rect() + 
#     scale_fill_viridis_c(na.value="grey85", limits=lims$uv_std) + 
#     facet_grid(.~lon) + 
#     labs(title=timesteps[i], x="Depth (m)", y="Latitude") + 
#     scale_x_continuous(breaks=c(0, 2000, 4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lat_slices_uv_std-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=6, units="in", dpi=300)
#   
# }
# 
# rm(latSlice_df); gc()



# * lon slices ------------------------------------------------------------

# lonSlice_df <- map_dfr(dir("data/HYCOM", "lon_slices_", full.names=T), readRDS) |>
#   group_by(depth) |>
#   mutate(u_std=c(scale(u, center=F)),
#          v_std=c(scale(v, center=F)),
#          uv_std=c(scale(uv, center=F))) |>
#   ungroup() 
# gc()
# 
# lims <- lonSlice_df |>
#   select(u, v, uv, u_std, v_std, uv_std) |>
#   reframe(across(everything(), ~range(.x, na.rm=T)))
# 
# timesteps <- unique(lonSlice_df$time)
# 
# for(i in seq_along(timesteps)) {
#   df_i <- lonSlice_df |> 
#     filter(time==timesteps[i]) |>
#     mutate(lon0=lon-lon_diff,
#            lon1=lon+lon_diff)
#   
#   # * * u -------------------------------------------------------------------
#   p <- df_i |>
#     ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth0, ymax=-depth1, fill=u)) + 
#     geom_rect() + 
#     scale_fill_continuous_divergingx(palette="RdBu", na.value="grey85",
#                                      mid=0, limits=lims$u) + 
#     facet_grid(latFac~.) + 
#     labs(title=timesteps[i], y="Depth (m)", x="Latitude")+ 
#     scale_y_continuous(breaks=c(0, -2000, -4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lon_slices_u-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=10, units="in", dpi=300)
#   
#   # * * v -------------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth0, ymax=-depth1, fill=v)) + 
#     geom_rect() + 
#     scale_fill_continuous_divergingx(palette="PRGn", na.value="grey85", 
#                                      mid=0, limits=lims$v) + 
#     facet_grid(latFac~.) + 
#     labs(title=timesteps[i], y="Depth (m)", x="Latitude") + 
#     scale_y_continuous(breaks=c(0, -2000, -4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lon_slices_v-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=10, units="in", dpi=300)
#   
#   # * * uv ------------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth0, ymax=-depth1, fill=uv)) + 
#     geom_rect() + 
#     scale_fill_viridis_c(na.value="grey85", limits=lims$uv) + 
#     facet_grid(latFac~.) + 
#     labs(title=timesteps[i], y="Depth (m)", x="Latitude") + 
#     scale_y_continuous(breaks=c(0, -2000, -4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lon_slices_uv-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=10, units="in", dpi=300)
#   
#   # * * uvDir ---------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth0, ymax=-depth1, fill=uvDir)) + 
#     geom_rect() + 
#     scale_fill_gradientn(colours=cmr$seasons, limits=c(-pi, pi),
#                          breaks=c(-pi, -pi/2, 0, pi/2, pi), 
#                          labels=c("W", "S", "E", "N", "W")) + 
#     facet_grid(latFac~.) + 
#     labs(title=timesteps[i], y="Depth (m)", x="Latitude") + 
#     scale_y_continuous(breaks=c(0, -2000, -4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lon_slices_uvDir-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=10, units="in", dpi=300)
#   
#   # * * u_std ---------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth0, ymax=-depth1, fill=u_std)) + 
#     geom_rect() + 
#     scale_fill_continuous_divergingx(palette="RdBu", na.value="grey85", 
#                                      mid=0, limits=lims$u_std) + 
#     facet_grid(latFac~.) + 
#     labs(title=timesteps[i], y="Depth (m)", x="Latitude") + 
#     scale_y_continuous(breaks=c(0, -2000, -4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lon_slices_u_std-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=10, units="in", dpi=300)
#   
#   
#   # * * v_std ---------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth0, ymax=-depth1, fill=v_std)) + 
#     geom_rect() + 
#     scale_fill_continuous_divergingx(palette="PRGn", na.value="grey85",
#                                      mid=0, limits=lims$v_std) + 
#     facet_grid(latFac~.) + 
#     labs(title=timesteps[i], y="Depth (m)", x="Latitude") + 
#     scale_y_continuous(breaks=c(0, -2000, -4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lon_slices_v_std-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=10, units="in", dpi=300)
#   
#   # * * uv_std --------------------------------------------------------------
#   p <- df_i |> 
#     ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth0, ymax=-depth1, fill=uv_std)) + 
#     geom_rect() + 
#     scale_fill_viridis_c(na.value="grey85", limits=lims$uv_std) + 
#     facet_grid(latFac~.) + 
#     labs(title=timesteps[i], y="Depth (m)", x="Latitude") + 
#     scale_y_continuous(breaks=c(0, -2000, -4000)) 
#   ggsave(glue("figs/anim/temp/HYCOM_lon_slices_uv_std-{str_pad(i,4,'left','0')}.png"), 
#          p, width=9, height=10, units="in", dpi=300)
#   
# }




# depth slices ------------------------------------------------------------


for(i in dirf(glue("{E_dir}data/HYCOM"), "dep_slices_")) {
  i_df <- readRDS(i)
  depths <- unique(i_df$depth)
  for(j in depths) {
    i_df |>
      filter(depth==j) |>
      saveRDS(glue("{str_replace(i, 'slices_', glue('slices_{j}_'))}"))
  }
  gc()
}

for(j in depths) {
  map_dfr(dirf("{E_dir}data/HYCOM", glue("dep_slices_{j}_")), readRDS) |>
    mutate(u_std=c(scale(u, center=F)),
           v_std=c(scale(v, center=F)),
           uv_std=c(scale(uv, center=F))) |>
    saveRDS(glue("data/HYCOM/depMerged_slices_{j}.rds"))
}
gc()

lims <- map_dfr(dirf("data/HYCOM", "depMerged_TEMP_slices"),
                  ~readRDS(.x) |>
                  select(depth, u, v, uv, u_std, v_std, uv_std) |>
                  reframe(across(matches("u|v"), ~range(.x, na.rm=T)), .by=depth))
gc()
lims <- lims |>
  reframe(across(matches("u|v"), ~range(.x)))
timesteps <- unique(readRDS(dirf("data/HYCOM", "depMerged_TEMP_slices")[1])$time)

theme_set(theme_bw() + theme(legend.position="bottom",
                             legend.key.width=unit(2, "cm"),
                             legend.key.height=unit(0.25, "cm")))
for(i in seq_along(timesteps)[-(1:8)]) {
  df_i <- map_dfr(dirf("data/HYCOM", "depMerged_TEMP_slices"),
                  ~readRDS(.x) |> filter(time==timesteps[i]))
  gc()

  # * * u -------------------------------------------------------------------
  # p <- df_i |>
  #   ggplot() +
  #   geom_raster(aes(lon-360, lat, fill=u)) +
  #   geom_sf(data=ccz, colour="grey30", linewidth=0.25) +
  #   scale_fill_continuous_divergingx(palette="RdBu", na.value="grey85",
  #                                    mid=0, limits=lims$u) +
  #   facet_grid(depth~.) +
  #   labs(title=timesteps[i], x="Longitude", y="Latitude")
  # ggsave(glue("figs/anim/temp/HYCOM_dep_slices_u-{str_pad(i,4,'left','0')}.png"),
  #        p, width=5, height=10, units="in", dpi=300)

  # * * v -------------------------------------------------------------------
  # p <- df_i |>
  #   ggplot() +
  #   geom_raster(aes(lon-360, lat, fill=v)) +
  #   geom_sf(data=ccz, colour="grey30", linewidth=0.25) +
  #   scale_fill_continuous_divergingx(palette="PRGn", na.value="grey85",
  #                                    mid=0, limits=lims$v) +
  #   facet_grid(depth~.) +
  #   labs(title=timesteps[i], x="Longitude", y="Latitude")
  # ggsave(glue("figs/anim/temp/HYCOM_dep_slices_v-{str_pad(i,4,'left','0')}.png"),
  #        p, width=5, height=10, units="in", dpi=300)

  # * * uv ------------------------------------------------------------------
  # p <- df_i |>
  #   ggplot() +
  #   geom_raster(aes(lon-360, lat, fill=uv)) +
  #   geom_sf(data=ccz, colour="grey30", linewidth=0.25) +
  #   scale_fill_viridis_c(na.value="grey85", limits=lims$uv) +
  #   facet_grid(depth~.) +
  #   labs(title=timesteps[i], x="Longitude", y="Latitude")
  # ggsave(glue("figs/anim/temp/HYCOM_dep_slices_uv-{str_pad(i,4,'left','0')}.png"),
  #        p, width=5, height=10, units="in", dpi=300)

  # * * uvDir ---------------------------------------------------------------
  p <- df_i |>
    ggplot() +
    geom_raster(aes(lon-360, lat, fill=uvDir)) +
    geom_sf(data=ccz, colour="grey30", linewidth=0.25) +
    scale_fill_gradientn(colours=cmr$seasons, limits=c(-pi, pi),
                         breaks=c(-pi, -pi/2, 0, pi/2, pi),
                         labels=c("W", "S", "E", "N", "W")) +
    facet_grid(depth~.) +
    labs(title=timesteps[i], x="Longitude", y="Latitude")
  ggsave(glue("figs/anim/temp/HYCOM_TEMP_dep_slices_uvDir-{str_pad(i,4,'left','0')}.png"),
         p, width=5, height=10, units="in", dpi=300)

  # * * u_std ---------------------------------------------------------------
  p <- df_i |>
    ggplot() +
    geom_raster(aes(lon-360, lat, fill=u_std)) +
    geom_sf(data=ccz, colour="grey30", linewidth=0.25) +
    scale_fill_continuous_divergingx(palette="RdBu", na.value="grey85",
                                     mid=0, limits=lims$u_std) +
    facet_grid(depth~.) +
    labs(title=timesteps[i], x="Longitude", y="Latitude")
  ggsave(glue("figs/anim/temp/HYCOM_TEMP_dep_slices_u_std-{str_pad(i,4,'left','0')}.png"),
         p, width=5, height=10, units="in", dpi=300)

  # * * v_std ---------------------------------------------------------------
  p <- df_i |>
    ggplot() +
    geom_raster(aes(lon-360, lat, fill=v_std)) +
    geom_sf(data=ccz, colour="grey30", linewidth=0.25) +
    scale_fill_continuous_divergingx(palette="PRGn", na.value="grey85",
                                     mid=0, limits=lims$v_std) +
    facet_grid(depth~.) +
    labs(title=timesteps[i], x="Longitude", y="Latitude")
  ggsave(glue("figs/anim/temp/HYCOM_TEMP_dep_slices_v_std-{str_pad(i,4,'left','0')}.png"),
         p, width=5, height=10, units="in", dpi=300)

  # * * uv_std --------------------------------------------------------------
  p <- df_i |>
    ggplot() +
    geom_raster(aes(lon-360, lat, fill=uv_std)) +
    geom_sf(data=ccz, colour="grey30", linewidth=0.25) +
    scale_fill_viridis_c(na.value="grey85", limits=lims$uv_std) +
    facet_grid(depth~.) +
    labs(title=timesteps[i], x="Longitude", y="Latitude")
  ggsave(glue("figs/anim/temp/HYCOM_TEMP_dep_slices_uv_std-{str_pad(i,4,'left','0')}.png"),
         p, width=5, height=10, units="in", dpi=300)

}

rm(i_df); gc()





# make gifs ---------------------------------------------------------------

# source("code/0000_make_hycom_gifs.R")


# deprecated --------------------------------------------------------------





# 
# 
# 
# anim <- lonSlice_df |> 
#   ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth1, ymax=-depth0, fill=u)) + 
#   geom_rect() + 
#   scale_fill_continuous_divergingx(palette="RdBu", na.value="grey85", mid=0) + 
#   facet_grid(latFac~.) + 
#   transition_manual(time) + 
#   labs(title="{current_frame}", x="Longitude", y="Depth (m)") + 
#   scale_y_continuous(breaks=c(0, -2000, -4000)) + 
#   theme(legend.position="bottom",
#         legend.key.width=unit(2, "cm"),
#         legend.key.height=unit(0.25, "cm"))
# anim_save("figs/anim/00_HYCOM_lon_slices_u.gif", anim, 
#           nframes=n_distinct(lonSlice_df$time), fps=25, 
#           width=9, height=10, units="in", res=300)
# 
# anim <- lonSlice_df |> 
#   ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth1, ymax=-depth0, fill=v)) + 
#   geom_rect() + 
#   scale_fill_continuous_divergingx(palette="PRGn", na.value="grey85", mid=0) + 
#   facet_grid(latFac~.) + 
#   transition_manual(time) + 
#   labs(title="{current_frame}", x="Longitude", y="Depth (m)") + 
#   scale_y_continuous(breaks=c(0, -2000, -4000)) + 
#   theme(legend.position="bottom",
#         legend.key.width=unit(2, "cm"),
#         legend.key.height=unit(0.25, "cm"))
# anim_save("figs/anim/00_HYCOM_lon_slices_v.gif", anim, 
#           nframes=n_distinct(lonSlice_df$time), fps=25, 
#           width=9, height=10, units="in", res=300)
# 
# anim <- lonSlice_df |> 
#   ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth1, ymax=-depth0, fill=uv)) + 
#   geom_rect() + 
#   scale_fill_viridis_c(na.value="grey85") + 
#   facet_grid(latFac~.) + 
#   transition_manual(time) + 
#   labs(title="{current_frame}", x="Longitude", y="Depth (m)") + 
#   scale_y_continuous(breaks=c(0, -2000, -4000)) + 
#   theme(legend.position="bottom",
#         legend.key.width=unit(2, "cm"),
#         legend.key.height=unit(0.25, "cm"))
# anim_save("figs/anim/00_HYCOM_lon_slices_uv.gif", anim, 
#           nframes=n_distinct(lonSlice_df$time), fps=25, 
#           width=9, height=10, units="in", res=300)
# 
# anim <- lonSlice_df |> 
#   ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth1, ymax=-depth0, fill=uvDir)) + 
#   geom_rect() + 
#   scale_fill_gradientn(colours=cmr$seasons, limits=c(-pi, pi),
#                        breaks=c(-pi, -pi/2, 0, pi/2, pi), 
#                        labels=c("W", "S", "E", "N", "W")) + 
#   facet_grid(latFac~.) + 
#   transition_manual(time) + 
#   labs(title="{current_frame}", x="Longitude", y="Depth (m)") + 
#   scale_y_continuous(breaks=c(0, -2000, -4000)) + 
#   theme(legend.position="bottom",
#         legend.key.width=unit(2, "cm"),
#         legend.key.height=unit(0.25, "cm"))
# anim_save("figs/anim/00_HYCOM_lon_slices_uvDir.gif", anim, 
#           nframes=n_distinct(lonSlice_df$time), fps=25, 
#           width=9, height=10, units="in", res=300)
# 
# 
# anim <- lonSlice_df |> 
#   ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth1, ymax=-depth0, fill=u_std)) + 
#   geom_rect() + 
#   scale_fill_continuous_divergingx(palette="RdBu", na.value="grey85", mid=0) + 
#   facet_grid(latFac~.) + 
#   transition_manual(time) + 
#   labs(title="{current_frame}", x="Longitude", y="Depth (m)") + 
#   scale_y_continuous(breaks=c(0, -2000, -4000)) + 
#   theme(legend.position="bottom",
#         legend.key.width=unit(2, "cm"),
#         legend.key.height=unit(0.25, "cm"))
# anim_save("figs/anim/00_HYCOM_lon_slices_u_std.gif", anim, 
#           nframes=n_distinct(lonSlice_df$time), fps=25, 
#           width=9, height=10, units="in", res=300)
# 
# anim <- lonSlice_df |> 
#   ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth1, ymax=-depth0, fill=v_std)) + 
#   geom_rect() + 
#   scale_fill_continuous_divergingx(palette="PRGn", na.value="grey85", mid=0) + 
#   facet_grid(latFac~.) + 
#   transition_manual(time) + 
#   labs(title="{current_frame}", x="Longitude", y="Depth (m)") + 
#   scale_y_continuous(breaks=c(0, -2000, -4000)) + 
#   theme(legend.position="bottom",
#         legend.key.width=unit(2, "cm"),
#         legend.key.height=unit(0.25, "cm"))
# anim_save("figs/anim/00_HYCOM_lon_slices_v_std.gif", anim, 
#           nframes=n_distinct(lonSlice_df$time), fps=25, 
#           width=9, height=10, units="in", res=300)
# 
# anim <- lonSlice_df |> 
#   ggplot(aes(xmin=lon0, xmax=lon1, ymin=-depth1, ymax=-depth0, fill=uv_std)) + 
#   geom_rect() + 
#   scale_fill_viridis_c(na.value="grey85") + 
#   facet_grid(latFac~.) + 
#   transition_manual(time) + 
#   labs(title="{current_frame}", x="Longitude", y="Depth (m)") + 
#   scale_y_continuous(breaks=c(0, -2000, -4000)) + 
#   theme(legend.position="bottom",
#         legend.key.width=unit(2, "cm"),
#         legend.key.height=unit(0.25, "cm"))
# anim_save("figs/anim/00_HYCOM_lon_slices_uv_std.gif", anim, 
#           nframes=n_distinct(lonSlice_df$time), fps=25, 
#           width=9, height=10, units="in", res=300)
# 
# 



