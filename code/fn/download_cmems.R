# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


download_cmems <- function(creds, 
                           bbox, 
                           dates=c(ymd("2023-01-01"), today()),
                           out_nc=NULL) {
  
  if(dates[1] < "2019-12-01") {
    url <- paste0("https://", creds$userid, ":", creds$pw, "@", 
                  "my.cmems-du.eu/thredds/dodsC/",
                  "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D")
  } else {
    url <- paste0("https://", creds$userid, ":", creds$pw, "@", 
                  "nrt.cmems-du.eu/thredds/dodsC/",
                  "dataset-duacs-nrt-global-merged-allsat-phy-l4") 
  }
  nc <- nc_open(url)
  
  # get metadata and identify lon/lat/time indexes to extract
  cmems_i <- list(
    longitude=ncvar_get(nc, "longitude"),
    latitude=ncvar_get(nc, "latitude"),
    time=ncvar_get(nc, "time")
  )
  cmems_i$date <- as_date(ymd(str_split_fixed(ncatt_get(nc, "time", "units")$value, " ", 4)[3]) + 
                            ddays(cmems_i$time))
  nc_close(nc)
  chunk_i <- list(
    longitude=between(cmems_i$longitude, bbox$xmin, bbox$xmax) |>
      which() |> range(),
    latitude=between(cmems_i$latitude, bbox$ymin, bbox$ymax) |>
      which() |> range(),
    time=between(cmems_i$date, dates[1], dates[2]) |> which() |> range()
  ) |>
    map(~paste0("[", min(.x)-1, ":", max(.x)-1, "]"))
  
  suffix <- paste0("?",
                   "longitude", chunk_i$longitude, ",",
                   "latitude", chunk_i$latitude, ",",
                   "time", chunk_i$time, ",",
                   "sla", chunk_i$time, chunk_i$latitude, chunk_i$longitude, ",",
                   "ugos", chunk_i$time, chunk_i$latitude, chunk_i$longitude, ",",
                   "vgos", chunk_i$time, chunk_i$latitude, chunk_i$longitude, ",",
                   "ugosa", chunk_i$time, chunk_i$latitude, chunk_i$longitude, ",",
                   "vgosa", chunk_i$time, chunk_i$latitude, chunk_i$longitude)
  
  # Generate URL and attempt to download
  nc_url <- paste0(url, suffix)
  cmems_nc <- nc_open(nc_url)
  cmems_names <- c(names(cmems_nc$dim), names(cmems_nc$var))
  cmems_ls <- map(cmems_names, 
                  ~ncvar_get(cmems_nc, .x)) |>
    set_names(cmems_names)
  
  # Get attributes
  atts <- cmems_names |>
    as.list() |>
    append(0) |>
    set_names(c(names(cmems_nc$dim), names(cmems_nc$var), "global")) |>
    map(~ncatt_get(cmems_nc, .x))
  # Define the dimensions
  dims <- list("longitude", "latitude", "time") |>
    set_names() |>
    map(~ncdim_def(.x, atts[[.x]]$units, cmems_ls[[.x]]))
  # Define variables
  vars <- names(cmems_nc$var) |>
    as.list() |> 
    set_names() |>
    map(~ncvar_def(.x, atts[[.x]]$units, dims))
  nc_close(cmems_nc) 
  
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
  walk(names(vars), ~ncvar_put(nc=ncnew, varid=vars[[.x]], vals=cmems_ls[[.x]]))
  for(v in names(vars)) {
    for(att in names(atts[[v]])) {
      ncatt_put(ncnew, v, att, atts[[v]][[att]])
    }
  }
  nc_close(ncnew) 
  
  
  cmems_df <- expand_grid(time=cmems_ls$time,
                          lat=cmems_ls$latitude,
                          lon=cmems_ls$longitude) |>
    mutate(date=as_date(ymd("1950-01-01") + ddays(time)),
           sla=c(cmems_ls$sla),
           ugos=c(cmems_ls$ugos),
           vgos=c(cmems_ls$vgos),
           ugosa=c(cmems_ls$ugosa),
           vgosa=c(cmems_ls$vgosa)) |>
    mutate(uv=sqrt(ugos^2 + vgos^2),
           EKE_cm2s2=((ugosa*100)^2 + (vgosa*100)^2) / 2)
  
  
  # nc.lon <- ncvar_get(nc, "longitude")
  # nc.lat <- ncvar_get(nc, "latitude")
  # nc.time <- ncvar_get(nc, "time")
  # time_units <- ncatt_get(nc, "time", "units")
  # nc.origin <- ymd(str_split_fixed(time_units$value, " ", 4)[3])
  # nc.date <- as_date(nc.origin + ddays(nc.time))
  
  # lon_i <- which(between(nc.lon, bbox$xmin, bbox$xmax))
  # lat_i <- which(between(nc.lat, bbox$ymin, bbox$ymax))
  # time_i <- which(between(nc.date, dates[1], dates[2]))
  # var.start <- c(lon_i[1], lat_i[1], time_i[1])
  # var.count <- c(length(lon_i), length(lat_i), length(time_i)) 
  # 
  # # extract variables
  # cmems_df <- expand_grid(date=nc.date[time_i],
  #                         lat=nc.lat[lat_i],
  #                         lon=nc.lon[lon_i]) |>
  #   mutate(sla=c(ncvar_get(nc, "sla", start=var.start, count=var.count)),
  #          ugos=c(ncvar_get(nc, "ugos", start=var.start, count=var.count)),
  #          vgos=c(ncvar_get(nc, "vgos", start=var.start, count=var.count)),
  #          ugosa=c(ncvar_get(nc, "ugosa", start=var.start, count=var.count)),
  #          vgosa=c(ncvar_get(nc, "vgosa", start=var.start, count=var.count))) |>
  #   mutate(uv=sqrt(ugos^2 + vgos^2),
  #          EKE_cm2s2=((ugosa*100)^2 + (vgosa*100)^2) / 2)
  
  # nc_close(nc)
  
  return(cmems_df)
  
}
