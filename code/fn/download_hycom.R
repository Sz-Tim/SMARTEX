# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


download_hycom <- function(hy_i, bbox, time_range, depth_range=c(1, 40), 
                           out_nc=NULL, out_rds=NULL, out_rng=NULL) {
  library(tidyverse); library(ncdf4)
  
  # Define subset indices
  chunk_i <- list(
    lon=between(hy_i$lon, 360+bbox$xmin, 360+bbox$xmax) |>
      which() |> range(),
    lat=between(hy_i$lat, bbox$ymin, bbox$ymax) |>
      which() |> range(),
    time=between(hy_i$time, time_range[1], time_range[2]) |> which() |> range(),
    depth=depth_range) |>
    map(~paste0("[", min(.x)-1, ":", max(.x)-1, "]"))
  chunk_suf <- paste0("?",
                      "depth", chunk_i$depth, ",",
                      "lon", chunk_i$lon, ",",
                      "lat", chunk_i$lat, ",",
                      "time", chunk_i$time, ",",
                      "water_u", chunk_i$time, chunk_i$depth, chunk_i$lat, chunk_i$lon, ",",
                      "water_v", chunk_i$time, chunk_i$depth, chunk_i$lat, chunk_i$lon)
  chunk_ls <- chunk_df.ls <- map(seq_along(chunk_suf), ~list())
  
  # Generate URL and attempt to download
  nc_url <- paste0(thredds_url, chunk_suf)
  hy_nc <- nc_open(nc_url)
  hy_ls <- list(
    depth=ncvar_get(hy_nc, "depth"),
    time=ncvar_get(hy_nc, "time"),
    lon=ncvar_get(hy_nc, "lon"),
    lat=ncvar_get(hy_nc, "lat"),
    u=ncvar_get(hy_nc, "water_u"),
    v=ncvar_get(hy_nc, "water_v")
  )
  
  if(!is.null(out_nc)) {
    # Define the dimensions
    dim_depth <- ncdim_def("depth", "m", hy_ls$depth)
    dim_time <- ncdim_def("time", "hours since 2000-01-01 00:00:00", hy_ls$time)
    dim_lat <- ncdim_def("lat", "degrees_north", hy_ls$lat)
    dim_lon <- ncdim_def("lon", "degrees_east", hy_ls$lon)
    
    # Define variables
    var_u <- ncvar_def("water_u", "m/s", 
                       list(dim_time, dim_depth, dim_lat, dim_lon), -30000, 
                       longname = "Eastward Water Velocity", prec = "float")
    var_v <- ncvar_def("water_v", "m/s", 
                       list(dim_time, dim_depth, dim_lat, dim_lon), -30000, 
                       longname = "Northward Water Velocity", prec = "float")
    
    # Create nc file
    ncnew <- nc_create(out_nc, list(var_u, var_v))
    ncvar_put(ncnew, var_u, hy_ls$u)
    ncvar_put(ncnew, var_v, hy_ls$v)
    nc_close(ncnew) 
  }
  
  if(!is.null(out_rds)) {
    # define midpoints between depth layers
    depth_layers <- tibble(depth=unique(hy_ls$depth)) |>
      mutate(depth0=depth - (depth-lag(depth))/2,
             depth1=depth + (lead(depth)-depth)/2) |>
      mutate(depth0=if_else(is.na(depth0), min(depth), depth0),
             depth1=if_else(is.na(depth1), max(depth), depth1))
    # create dataframe
    hy_df <- expand_grid(time=hy_ls$time,
                               depth=hy_ls$depth,
                               lat=hy_ls$lat,
                               lon=hy_ls$lon) |>
      arrange(time, depth, lat, lon) |>
      mutate(id=as.numeric(as.factor(paste(lon, lat))),
             time=as_datetime("2000-01-01 00:00:00") + time*60*60,
             u=c(hy_ls$u),
             v=c(hy_ls$v)) |>
      mutate(u2v2=u^2 + v^2,
             uv=sqrt(u2v2),
             uvDir=atan2(v,u)) |>
      left_join(depth_layers)
    saveRDS(hy_df, out_rds)
  }
  if(!is.null(out_rng)) {
    if(is.null(out_rds)) { 
      return("Must provide out_rds for range calculation") 
    }
    ranges <- hy_df |>
      group_by(depth) |>
      summarise(u_min=min(u, na.rm=T),
                u_max=max(u, na.rm=T),
                v_min=min(v, na.rm=T),
                v_max=max(v, na.rm=T),
                u2v2_min=min(u2v2, na.rm=T),
                u2v2_max=max(u2v2, na.rm=T),
                uv_min=min(uv, na.rm=T),
                uv_max=max(uv, na.rm=T),
                u_mean=mean(u, na.rm=T),
                uv_mean=mean(v, na.rm=T),
                u2v2_mean=mean(u2v2, na.rm=T),
                uv_mean=mean(uv, na.rm=T),
                u_sd=sd(u, na.rm=T),
                v_sd=sd(v, na.rm=T),
                u2v2_sd=sd(u2v2, na.rm=T),
                uv_sd=sd(uv, na.rm=T)) |>
      ungroup()
    saveRDS(ranges, out_rng)
  }
  return(cat("Created: \n  ", out_nc, "\n  ", out_rds, "\n  ", out_rng, "\n"))
}
