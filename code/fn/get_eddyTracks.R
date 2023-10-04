# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk



get_eddyTracks <- function(bbox, dates=c(ymd("2023-01-01"), today())) {
  
  track_df <- vector("list", 2) |> setNames(c("anticyclonic", "cyclonic"))
  for(i in names(track_df)) {
    nc <- nc_open(dirf("data/AVISO", glue("_{i}.*nc")))
    nc.lon <- ncvar_get(nc, "longitude")
    nc.lon[nc.lon > 180] <- nc.lon[nc.lon > 180] - 360
    nc.lat <- ncvar_get(nc, "latitude")
    nc.time <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")
    nc.origin <- ymd(str_split_fixed(time_units$value, " ", 4)[3])
    nc.date <- as_date(nc.origin + ddays(nc.time))
    lon_i <- between(nc.lon, bbox$xmin, bbox$xmax)
    lat_i <- between(nc.lat, bbox$ymin, bbox$ymax)
    time_i <- between(nc.date, dates[1], dates[2])
    rows_inbound <- which(lon_i & lat_i & time_i)
    
    track_df[[i]] <- tibble(
      type=i,
      aviso_id=rows_inbound,
      lon=nc.lon[rows_inbound],
      lat=nc.lat[rows_inbound],
      date=nc.date[rows_inbound],
      amplitude=ncvar_get(nc, "amplitude")[rows_inbound],
      effective_radius=ncvar_get(nc, "effective_radius")[rows_inbound],
      track=ncvar_get(nc, "track")[rows_inbound],
      observation_number=ncvar_get(nc, "observation_number")[rows_inbound]
    ) |>
      arrange(track, observation_number)
    nc_close(nc)
  }
  track_df <- reduce(track_df, bind_rows)
  return(track_df)
}