# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


download_hycom <- function(url_base="http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0", 
                           hy_i, 
                           bbox, time_range, depth_range, 
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
                      "water_v", chunk_i$time, chunk_i$depth, chunk_i$lat, chunk_i$lon, ",",
                      "water_temp", chunk_i$time, chunk_i$depth, chunk_i$lat, chunk_i$lon, ",",
                      "salinity", chunk_i$time, chunk_i$depth, chunk_i$lat, chunk_i$lon, ",",
                      "surf_el", chunk_i$time, chunk_i$lat, chunk_i$lon, ",",
                      "water_u_bottom", chunk_i$time, chunk_i$lat, chunk_i$lon, ",",
                      "water_v_bottom", chunk_i$time, chunk_i$lat, chunk_i$lon, ",",
                      "water_temp_bottom", chunk_i$time, chunk_i$lat, chunk_i$lon, ",",
                      "salinity_bottom", chunk_i$time, chunk_i$lat, chunk_i$lon)
  
  # Generate URL and attempt to download
  nc_url <- paste0(url_base, chunk_suf)
  if(file.exists(out_nc)) {
    nc_url <- out_nc
  }
  hy_nc <- nc_open(nc_url)
  hy_ls <- list(
    time=ncvar_get(hy_nc, "time"),
    lon=ncvar_get(hy_nc, "lon"),
    lat=ncvar_get(hy_nc, "lat")
  )
  hy_ls$depth <- ncvar_get(hy_nc, "depth")
  hy_ls$water_u <- ncvar_get(hy_nc, "water_u")
  hy_ls$water_v <- ncvar_get(hy_nc, "water_v")
  hy_ls$water_temp <- ncvar_get(hy_nc, "water_temp")
  hy_ls$salinity <- ncvar_get(hy_nc, "salinity")
  hy_ls$surf_el <- ncvar_get(hy_nc, "surf_el")
  hy_ls$water_u_bottom <- ncvar_get(hy_nc, "water_u_bottom")
  hy_ls$water_v_bottom <- ncvar_get(hy_nc, "water_v_bottom")
  hy_ls$water_temp_bottom <- ncvar_get(hy_nc, "water_temp_bottom")
  hy_ls$salinity_bottom <- ncvar_get(hy_nc, "salinity_bottom")
  
  # Get attributes
  atts <- c(names(hy_nc$dim), names(hy_nc$var)) |>
    as.list() |>
    append(0) |>
    set_names(c(names(hy_nc$dim), names(hy_nc$var), "global")) |>
    map(~ncatt_get(hy_nc, .x))
  # Define the dimensions
  dims <- list("lon", "lat", "depth", "time") |>
    discard(is.na) |>
    set_names() |>
    map(~ncdim_def(.x, atts[[.x]]$units, hy_ls[[.x]]))
  # Define variables
  dims_ls <- list(
    xyz=dims,
    xy=dims[-3]
  )
  vars <- names(hy_nc$var) |>
    as.list() |> 
    set_names() |>
    map(~ncvar_def(.x, atts[[.x]]$units, dims_ls[[(grepl("surf|bottom", .x))+1]]))
  nc_close(hy_nc)
  
  # Store nc as out_nc
  if(!is.null(out_nc) & !file.exists(out_nc)) {
    # Create nc file
    ncnew <- nc_create(out_nc, vars)
    # global attributes
    for(att in names(atts[["global"]])) {
      ncatt_put(ncnew, 0, att, atts[["global"]][[att]])
    }
    # dimension attributes
    for(d in names(dims)) {
      for(att in names(atts[[d]])) {
        ncatt_put(ncnew, d, att, atts[[d]][[att]])
      }
    }
    # variables
    walk(names(vars), ~ncvar_put(nc=ncnew, varid=vars[[.x]], vals=hy_ls[[.x]]))
    for(v in names(vars)) {
      for(att in names(atts[[v]])) {
        ncatt_put(ncnew, v, att, atts[[v]][[att]])
      }
    }
    nc_close(ncnew) 
  }
  
  if(!is.null(out_rds)) {
      # create dataframe
      hy_df <- expand_grid(time=hy_ls$time,
                           depth=hy_ls$depth,
                           lat=hy_ls$lat,
                           lon=hy_ls$lon) |>
        arrange(time, depth, lat, lon) |>
        mutate(id=as.numeric(as.factor(paste(lon, lat))),
               time=as_datetime("2000-01-01 00:00:00") + time*60*60,
               u=c(hy_ls$water_u),
               v=c(hy_ls$water_v),
               temperature=c(hy_ls$water_temp),
               salinity=c(hy_ls$salinity)) |>
        mutate(u2v2=u^2 + v^2,
               uv=sqrt(u2v2),
               uvDir=atan2(v,u))
      saveRDS(hy_df, str_replace(out_rds, ".rds", "_3D.rds"))
      hy_2D_df <- expand_grid(time=hy_ls$time,
                              lat=hy_ls$lat,
                              lon=hy_ls$lon) |>
        arrange(time, lat, lon) |>
        mutate(id=as.numeric(as.factor(paste(lon, lat))),
               time=as_datetime("2000-01-01 00:00:00") + time*60*60,
               surf_el=c(hy_ls$surf_el),
               water_u_bottom=c(hy_ls$water_u_bottom),
               water_v_bottom=c(hy_ls$water_v_bottom),
               water_temp_bottom=c(hy_ls$water_temp_bottom),
               salinity_bottom=c(hy_ls$salinity_bottom)) |>
        mutate(u2v2_bottom=water_u_bottom^2 + water_v_bottom^2,
               uv_bottom=sqrt(u2v2_bottom),
               uvDir_bottom=atan2(water_v_bottom, water_u_bottom))
    saveRDS(hy_2D_df, str_replace(out_rds, ".rds", "_2D.rds"))
  }
  if(!is.null(out_rng)) {
    ranges <- hy_df |>
      group_by(depth) |>
      summarise(u_min=min(u, na.rm=T),
                u_max=max(u, na.rm=T),
                u_mean=mean(u, na.rm=T),
                u_sd=sd(u, na.rm=T),
                
                v_min=min(v, na.rm=T),
                v_max=max(v, na.rm=T),
                v_mean=mean(v, na.rm=T),
                v_sd=sd(v, na.rm=T),
                
                u2v2_min=min(u2v2, na.rm=T),
                u2v2_max=max(u2v2, na.rm=T),
                u2v2_mean=mean(u2v2, na.rm=T),
                u2v2_sd=sd(u2v2, na.rm=T),
                
                uv_min=min(uv, na.rm=T),
                uv_max=max(uv, na.rm=T),
                uv_mean=mean(uv, na.rm=T),
                uv_sd=sd(uv, na.rm=T),
                
                temperature_min=min(temperature, na.rm=T),
                temperature_max=max(temperature, na.rm=T),
                temperature_mean=mean(temperature, na.rm=T),
                temperature_sd=sd(temperature, na.rm=T),
                
                salinity_min=min(salinity, na.rm=T),
                salinity_max=max(salinity, na.rm=T),
                salinity_mean=mean(salinity, na.rm=T),
                salinity_sd=sd(salinity, na.rm=T)) |>
      ungroup()
    saveRDS(ranges, str_replace(out_rng, ".rds", "_3D.rds"))
    ranges_2D <- hy_2D_df |>
      summarise(surf_el_min=min(surf_el, na.rm=T),
                surf_el_max=max(surf_el, na.rm=T),
                surf_el_mean=mean(surf_el, na.rm=T),
                surf_el_sd=sd(surf_el, na.rm=T),
                
                water_u_bottom_min=min(water_u_bottom, na.rm=T),
                water_u_bottom_max=max(water_u_bottom, na.rm=T),
                water_u_bottom_mean=mean(water_u_bottom, na.rm=T),
                water_u_bottom_sd=sd(water_u_bottom, na.rm=T),
                
                water_v_bottom_min=min(water_v_bottom, na.rm=T),
                water_v_bottom_max=max(water_v_bottom, na.rm=T),
                water_v_bottom_mean=mean(water_v_bottom, na.rm=T),
                water_v_bottom_sd=sd(water_v_bottom, na.rm=T),
                
                u2v2_bottom_min=min(u2v2_bottom, na.rm=T),
                u2v2_bottom_max=max(u2v2_bottom, na.rm=T),
                u2v2_bottom_mean=mean(u2v2_bottom, na.rm=T),
                u2v2_bottom_sd=sd(u2v2_bottom, na.rm=T),
                
                uv_bottom_min=min(uv_bottom, na.rm=T),
                uv_bottom_max=max(uv_bottom, na.rm=T),
                uv_bottom_mean=mean(uv_bottom, na.rm=T),
                uv_bottom_sd=sd(uv_bottom, na.rm=T),
                
                water_temp_bottom_min=min(water_temp_bottom, na.rm=T),
                water_temp_bottom_max=max(water_temp_bottom, na.rm=T),
                water_temp_bottom_mean=mean(water_temp_bottom, na.rm=T),
                water_temp_bottom_sd=sd(water_temp_bottom, na.rm=T),
                
                salinity_bottom_min=min(salinity_bottom, na.rm=T),
                salinity_bottom_max=max(salinity_bottom, na.rm=T),
                salinity_bottom_mean=mean(salinity_bottom, na.rm=T),
                salinity_bottom_sd=sd(salinity_bottom, na.rm=T)) |>
      ungroup()
    saveRDS(ranges_2D, str_replace(out_rng, ".rds", "_2D.rds"))
    
  }
  return(cat("Created: \n  ", out_nc, "\n  ", out_rds, "\n  ", out_rng, "\n"))
}





download_hycom_ftp <- function(url_base="http://data.hycom.org/datasets/GLBy0.08/expt_93.0/data/hindcasts", 
                               date_range, dir_nc) {
  library(xml2); library(rvest)
  date_seq <- seq(date_range[1], date_range[2], by=1)
  dateChr_seq <- str_remove_all(date_seq, "-")
  years <- unique(year(date_seq))
  
  for(i in seq_along(years)) {
    url_i <- glue("{url_base}/{years[i]}/")
    links <- read_html(url_i) |> 
      html_nodes("a") |>
      html_attr("href") |> 
      grep("ssh.nc", x=_, value=T) |>
      unique()
    link_dates <- str_sub(links, 16, 23)
    links_to_get <- links[link_dates %in% dateChr_seq]
    for(j in seq_along(links_to_get)) {
      if(j %% 20 == 0) {
        cat("Time: ", format(Sys.time()))
      }
      download.file(glue("{url_i}{links_to_get[j]}"), 
                    glue("{dir_nc}{links_to_get[j]}"), mode="wb")
    }
  }
  
}
