# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


download_cmems <- function(creds, 
                           bbox, 
                           dates=c(ymd("2023-01-01"), today())) {
  
  url <- paste0("https://", creds$userid, ":", creds$pw, "@", 
                "nrt.cmems-du.eu/thredds/dodsC/",
                "dataset-duacs-nrt-global-merged-allsat-phy-l4")
  nc <- nc_open(url)
  
  # get metadata and identify lon/lat/time indexes to extract
  nc.lon <- ncvar_get(nc, "longitude")
  nc.lat <- ncvar_get(nc, "latitude")
  nc.time <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")
  nc.origin <- ymd(str_split_fixed(time_units$value, " ", 4)[3])
  nc.date <- as_date(nc.origin + ddays(nc.time))
  lon_i <- which(between(nc.lon, bbox$xmin, bbox$xmax))
  lat_i <- which(between(nc.lat, bbox$ymin, bbox$ymax))
  time_i <- which(between(nc.date, dates[1], dates[2]))
  var.start <- c(lon_i[1], lat_i[1], time_i[1])
  var.count <- c(length(lon_i), length(lat_i), length(time_i)) 
  
  # extract variables
  cmems_df <- expand_grid(date=nc.date[time_i],
                          lat=nc.lat[lat_i],
                          lon=nc.lon[lon_i]) |>
    mutate(sla=c(ncvar_get(nc, "sla", start=var.start, count=var.count)),
           adt=c(ncvar_get(nc, "adt", start=var.start, count=var.count)),
           ugosa=c(ncvar_get(nc, "ugosa", start=var.start, count=var.count)),
           vgosa=c(ncvar_get(nc, "vgosa", start=var.start, count=var.count))) |>
    mutate(EKE_cm2s2=((ugosa*100)^2 + (vgosa*100)^2) / 2)
  
  return(cmems_df)
  
}
