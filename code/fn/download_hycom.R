# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


download_hycom <- function(url_base="http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z", 
                           hy_i, 
                           bbox, time_range, depth_range=c(1, 40), 
                           out_nc=NULL, out_rds=NULL, out_rng=NULL) {
  library(tidyverse); library(ncdf4)
  
  uv3z <- grepl("uv3z", url_base)
  
  # Define subset indices
  chunk_i <- list(
    lon=between(hy_i$lon, 360+bbox$xmin, 360+bbox$xmax) |>
      which() |> range(),
    lat=between(hy_i$lat, bbox$ymin, bbox$ymax) |>
      which() |> range(),
    time=between(hy_i$time, time_range[1], time_range[2]) |> which() |> range(),
    depth=depth_range) |>
    map(~paste0("[", min(.x)-1, ":", max(.x)-1, "]"))
  if(uv3z) {
    chunk_suf <- paste0("?",
                        "depth", chunk_i$depth, ",",
                        "lon", chunk_i$lon, ",",
                        "lat", chunk_i$lat, ",",
                        "time", chunk_i$time, ",",
                        "water_u", chunk_i$time, chunk_i$depth, chunk_i$lat, chunk_i$lon, ",",
                        "water_v", chunk_i$time, chunk_i$depth, chunk_i$lat, chunk_i$lon)
  } else {
    chunk_suf <- paste0("?",
                        "lon", chunk_i$lon, ",",
                        "lat", chunk_i$lat, ",",
                        "time", chunk_i$time, ",",
                        "surf_el", chunk_i$time, chunk_i$lat, chunk_i$lon)
  }
  
  
  # Generate URL and attempt to download
  nc_url <- paste0(url_base, chunk_suf)
  hy_nc <- nc_open(nc_url)
  hy_ls <- list(
    time=ncvar_get(hy_nc, "time"),
    lon=ncvar_get(hy_nc, "lon"),
    lat=ncvar_get(hy_nc, "lat")
  )
  if(uv3z) {
    hy_ls$depth <- ncvar_get(hy_nc, "depth")
    hy_ls$water_u <- ncvar_get(hy_nc, "water_u")
    hy_ls$water_v <- ncvar_get(hy_nc, "water_v")
  } else {
    hy_ls$surf_el <- ncvar_get(hy_nc, "surf_el")
  }
  
  # Get attributes
  atts <- c(names(hy_nc$dim), names(hy_nc$var)) |>
    as.list() |>
    append(0) |>
    set_names(c(names(hy_nc$dim), names(hy_nc$var), "global")) |>
    map(~ncatt_get(hy_nc, .x))
  # Define the dimensions
  dims <- list("lon", "lat", ifelse(uv3z, "depth", NA), "time") |>
    discard(is.na) |>
    set_names() |>
    map(~ncdim_def(.x, atts[[.x]]$units, hy_ls[[.x]]))
  # Define variables
  vars <- names(hy_nc$var) |>
    as.list() |> 
    set_names() |>
    map(~ncvar_def(.x, atts[[.x]]$units, dims))
  nc_close(hy_nc)
  
  if(!is.null(out_nc)) {
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
    walk(names(vars), ~ncvar_put(ncnew, vars[[.x]], hy_ls[[.x]]))
    for(v in names(vars)) {
      for(att in names(atts[[v]])) {
        ncatt_put(ncnew, v, att, atts[[v]][[att]])
      }
    }
    nc_close(ncnew) 
  }
  
  if(!is.null(out_rds)) {
    if(uv3z) {
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
    } else {
      hy_df <- expand_grid(time=hy_ls$time,
                           lat=hy_ls$lat,
                           lon=hy_ls$lon) |>
        arrange(time, lat, lon) |>
        mutate(id=as.numeric(as.factor(paste(lon, lat))),
               time=as_datetime("2000-01-01 00:00:00") + time*60*60,
               surf_el=c(hy_ls$surf_el))
    }
    saveRDS(hy_df, out_rds)
  }
  if(!is.null(out_rng)) {
    if(is.null(out_rds)) { 
      return("Must provide out_rds for range calculation") 
    }
    if(uv3z) {
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
    } else {
      ranges <- hy_df |>
        summarise(surf_el_min=min(surf_el, na.rm=T),
                  surf_el_max=max(surf_el, na.rm=T),
                  surf_el_mean=mean(surf_el, na.rm=T),
                  surf_el_sd=sd(surf_el, na.rm=T)) |>
        ungroup()
    }
    
    saveRDS(ranges, out_rng)
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
