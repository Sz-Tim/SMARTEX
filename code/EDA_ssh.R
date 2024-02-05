# Sea surface height
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
reload <- TRUE
POI <- read_csv("data/SMARTEX_locations_copy.csv") |>
  filter(grepl("BGR|Long|Atoll|001|002", Name_Full))

ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
east_bbox <- list(xmin=-130, xmax=-90, ymin=2, ymax=20)
ccz_bbox <- list(xmin=-160, xmax=-90, ymin=0, ymax=23.5)
domain_bbox <- list(xmin=-119, xmax=-114, ymin=12, ymax=15)
bathy <- read_stars("data/bathymetry/gebco_2023_n23.5_s0.0_w-160.0_e-90.0.tif")
land <- st_read("data/CCZ_areas/CCZ_box_coastline.gpkg")

dir_figTemp <- "D:/SMARTEX/figs/temp/"




# Images: Eastern focus 1993-2022 -----------------------------------------

if(reload) {
  for(i in 1993:2022) {
    cat("starting", i, "\n")
    cmems_yr_i <- download_cmems(creds=readRDS("data/cmems_cred.rds"),
                                 east_bbox, 
                                 dates=c(ymd(glue("{i}-01-01")), 
                                         ymd(glue("{i}-12-31"))),
                                 out_nc=glue("data/CMEMS/cmems_east_{i}.nc")) |>
      mutate(EKE=pmin(EKE_cm2s2, 3e3))
    saveRDS(cmems_yr_i, glue("data/CMEMS/cmems_east_{i}.rds"))
  }
}
cmems_east <- dirf("data/CMEMS", "cmems_east_[1-2][0-9][0-9][0-9].rds") |>
  map_dfr(readRDS)
gc()

# EKE 2D maps
EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
                out_f=c("mean", "sd", "variance"),
                legend=c("EKE mean", "EKE sd", "EKE variance"))[1,]
for(k in 1:nrow(EKE_i)) {
  plot_EKE_2D(cmems_east |> filter(!is.na(EKE)), 
              bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
              east_bbox[1:2], east_bbox[3:4],
              glue("figs/anim/east_EKE_{EKE_i$out_f[k]}_1993-2022.png"), 
              c(10, 5.5))
}
EKE_sf <- cmems_east |>
  group_by(lat, lon) |>
  summarise(EKE_mn=mean(EKE),
            EKE_var=var(EKE),
            EKE_sd=sd(EKE)) |>
  ungroup() |>
  st_as_sf(coords=c("lon", "lat"), crs=4326)
walk(EKE_i$var, ~EKE_sf |> select(all_of(.x), geometry) |>
       stars::st_rasterize() |>
       stars::write_stars(glue("data/CMEMS/{.x}_east_1993-2022.tif")))




# Images: CCZ 1993-2022 ----------------------------------------------

if(reload) {
  for(i in 1993:2022) {
    cat("starting", i, "\n")
    cmems_yr_i <- download_cmems(creds=readRDS("data/cmems_cred.rds"),
                                 ccz_bbox, 
                                 dates=c(ymd(glue("{i}-01-01")), 
                                         ymd(glue("{i}-12-31"))),
                                 out_nc=glue("data/CMEMS/cmems_ccz_{i}.nc")) |>
      mutate(EKE=pmin(EKE_cm2s2, 3e3))
    gc()
    saveRDS(cmems_yr_i, glue("data/CMEMS/cmems_ccz_{i}.rds"))
    gc()
  }
}
cmems_ccz <- dirf("data/CMEMS", "cmems_ccz_[1-2][0-9][0-9][0-9].rds") |>
  map_dfr(~readRDS(.x) |> filter(lat >= 3))
gc()

# EKE 2D maps
EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
                out_f=c("mean", "sd", "variance"),
                legend=c("EKE mean", "EKE sd", "EKE variance"))[1,]
for(k in 1:nrow(EKE_i)) {
  plot_EKE_2D(cmems_ccz |> filter(!is.na(EKE)), 
              bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
              ccz_bbox[1:2], ccz_bbox[3:4],
              glue("figs/anim/CCZ_EKE_{EKE_i$out_f[k]}_3N-23.5N_1993-2022.png"), 
              c(10, 5.5))
}
EKE_sf <- cmems_ccz |>
  group_by(lat, lon) |>
  summarise(EKE_mn=mean(EKE),
            EKE_var=var(EKE),
            EKE_sd=sd(EKE)) |>
  ungroup() |>
  st_as_sf(coords=c("lon", "lat"), crs=4326)
walk(EKE_i$var, ~EKE_sf |> select(all_of(.x), geometry) |>
       stars::st_rasterize() |>
       stars::write_stars(glue("data/CMEMS/{.x}_3N-23.5N_CCZ_1993-2022.tif")))


ccz_f <- dirf("data/CMEMS", "cmems_ccz_[1-2][0-9][0-9][0-9].rds")
ccz_rng <- ccz_f |>
  map_dfr(~readRDS(.x) |> 
            mutate(sqrt_abs_sla=sqrt(abs(sla))) |>
            reframe(across(5:13, ~range(.x, na.rm=T)))) |>
  reframe(across(everything(), ~range(.x, na.rm=T)))
cmems_vars <- list("sla"=ccz_rng$sla,
                   "EKE"=ccz_rng$EKE,
                   "uv"=ccz_rng$uv)
                   # "uv"=c(0,1))
annual_ls <- monthly_ls <- weekly_ls <- vector("list", length(ccz_f))
for(f in seq_along(ccz_f)) {
  cmems_f <- readRDS(ccz_f[f]) |>
    filter(lat > 3) |>
    mutate(year=year(date),
           month=month(date),
           week=week(date))
  
  annual_ls[[f]] <- cmems_f |>
    summarise(EKE=mean(EKE_cm2s2),
              uv=mean(uv),
              .by=c("lat", "lon", "year"))
  monthly_ls[[f]] <- cmems_f |>
    summarise(EKE=mean(EKE_cm2s2),
              uv=mean(uv),
              .by=c("lat", "lon", "year", "month"))
  weekly_ls[[f]] <- cmems_f |>
    summarise(EKE=mean(EKE_cm2s2),
              uv=mean(uv),
              .by=c("lat", "lon", "year", "week"))
  gc()
  cat("Finished", f, "\n")
}
annual_df <- reduce(annual_ls, bind_rows)
monthly_df <- reduce(monthly_ls, bind_rows) |>
  summarise(EKE_mn=mean(EKE), EKE_sd=sd(EKE),
            uv_mn=mean(uv), uv_sd=sd(uv),
            .by=c("lat", "lon", "month")) |>
  mutate(EKE_cv=EKE_sd/EKE_mn,
         uv_cv=uv_sd/uv_mn)
weekly_df <- reduce(weekly_ls, bind_rows) |>
  summarise(EKE_mn=mean(EKE), EKE_sd=sd(EKE), 
            uv_mn=mean(uv), uv_sd=sd(uv),
            .by=c("lat", "lon", "week")) |>
  mutate(EKE_cv=EKE_sd/EKE_mn,
         uv_cv=uv_sd/uv_mn)
rm(annual_ls); rm(monthly_ls); rm(weekly_ls); gc()

annual_rng <- annual_df |>
  reframe(across(4:5, ~range(.x, na.rm=T)))
for(y in unique(annual_df$year)) {
  plot_cmems(df=annual_df |> filter(year==y) |> filter(!is.na(EKE)),
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE", fill_lim=annual_rng$EKE,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=y, out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_mean_{y}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=annual_df |> filter(year==y) |> filter(!is.na(uv)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv", fill_lim=annual_rng$uv,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=y, out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_mean_year_{y}.png"),
             tracks=NULL, darkLines=F)
}

monthly_rng <- monthly_df |>
  reframe(across(4:9, ~range(.x, na.rm=T)))
for(m in 1:12) {
  plot_cmems(df=monthly_df |> filter(month==m) |> filter(!is.na(EKE_mn)),
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_mn", fill_lim=monthly_rng$EKE_mn,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=month.name[m], out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_mean_month-{str_pad(m, 2, 'left', '0')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=monthly_df |> filter(month==m) |> filter(!is.na(EKE_mn)),
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_sd", fill_lim=monthly_rng$EKE_sd,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=month.name[m], out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_sdInterannual_month-{str_pad(m, 2, 'left', '0')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=monthly_df |> filter(month==m) |> filter(!is.na(EKE_cv)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_cv", fill_lim=monthly_rng$EKE_cv,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=month.name[m], out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_cv_month-{str_pad(m, 2, 'left', '0')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=monthly_df |> filter(month==m) |> filter(!is.na(uv_mn)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_mn", fill_lim=monthly_rng$uv_mn,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=month.name[m], out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_mean_month-{str_pad(m, 2, 'left', '0')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=monthly_df |> filter(month==m) |> filter(!is.na(uv_sd)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_sd", fill_lim=monthly_rng$uv_sd,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=month.name[m], out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_sdInterannual_month-{str_pad(m, 2, 'left', '0')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=monthly_df |> filter(month==m) |> filter(!is.na(uv_cv)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_cv", fill_lim=monthly_rng$uv_cv,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=month.name[m], out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_cv_month-{str_pad(m, 2, 'left', '0')}.png"),
             tracks=NULL, darkLines=F)
}

weekly_rng <- weekly_df |>
  mutate(EKE_cv=EKE_sd/EKE_mn,
         uv_cv=uv_sd/uv_mn) |>
  reframe(across(4:9, ~range(.x, na.rm=T)))
for(w in 1:52) {
  date_w <- ymd("2023-01-01")+7*(w-1) 
  plot_cmems(df=weekly_df |> filter(week==w) |> filter(!is.na(EKE_mn)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_mn", fill_lim=weekly_rng$EKE_mn,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("Week of", format(date_w, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_mean_week-{format(date_w, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=weekly_df |> filter(week==w) |> filter(!is.na(EKE_sd)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_sd", fill_lim=weekly_rng$EKE_sd,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("Week of", format(date_w, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_sdInterannual_week-{format(date_w, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=weekly_df |> filter(week==w) |> filter(!is.na(EKE_cv)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_cv", fill_lim=weekly_rng$EKE_cv,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("Week of", format(date_w, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_cv_week-{format(date_w, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  
  plot_cmems(df=weekly_df |> filter(week==w) |> filter(!is.na(uv_mn)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_mn", fill_lim=weekly_rng$uv_mn,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("Week of", format(date_w, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_mean_week-{format(date_w, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=weekly_df |> filter(week==w) |> filter(!is.na(uv_sd)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_sd", fill_lim=weekly_rng$uv_sd,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("Week of", format(date_w, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_sdInterannual_week-{format(date_w, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=weekly_df |> filter(week==w) |> filter(!is.na(uv_cv)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_cv", fill_lim=weekly_rng$uv_cv,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("Week of", format(date_w, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_cv_week-{format(date_w, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
}


daily_df <- dirf("data/CMEMS", "cmems_ccz_[1-2][0-9][0-9][0-9].rds") |>
  map_dfr(~readRDS(.x) |> 
            filter(lat >= 3) |> 
            select(date, lat, lon, uv, EKE) |>
            mutate(yday=yday(date)) |> select(-date)) |>
  summarise(EKE_mn=mean(EKE), EKE_sd=sd(EKE), 
            uv_mn=mean(uv), uv_sd=sd(uv),
            .by=c("lat", "lon", "yday")) |>
  filter(yday < 366) |> 
  mutate(EKE_cv=EKE_sd/EKE_mn,
         uv_cv=uv_sd/uv_mn,
         date=ymd("2021-01-01") + yday-1) |> 
  select(-yday) 
saveRDS(daily_df, "data/CMEMS/CCZ_climatology_daily.rds")
daily_rng <- daily_df |>
  select(date, lon, lat, starts_with("uv"), starts_with("EKE")) |>
  reframe(across(4:9, ~range(.x, na.rm=T)))
for(y in 1:365) {
  date_y <- ymd("2021-01-01")+y-1 
  plot_cmems(df=daily_df |> filter(date==date_y) |> filter(!is.na(EKE_mn)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_mn", fill_lim=daily_rng$EKE_mn,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("", format(date_y, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_mean_day-{format(date_y, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=daily_df |> filter(date==date_y) |> filter(!is.na(EKE_sd)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_sd", fill_lim=daily_rng$EKE_sd,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("", format(date_y, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_sdInterannual_day-{format(date_y, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=daily_df |> filter(date==date_y) |> filter(!is.na(EKE_cv)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="EKE_cv", fill_lim=daily_rng$EKE_cv,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("", format(date_y, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_EKE_cv_day-{format(date_y, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  
  plot_cmems(df=daily_df |> filter(date==date_y) |> filter(!is.na(uv_mn)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_mn", fill_lim=daily_rng$uv_mn,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("", format(date_y, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_mean_day-{format(date_y, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=daily_df |> filter(date==date_y) |> filter(!is.na(uv_sd)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_sd", fill_lim=daily_rng$uv_sd,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("", format(date_y, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_sdInterannual_day-{format(date_y, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
  plot_cmems(df=daily_df |> filter(date==date_y) |> filter(!is.na(uv_cv)), 
             bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="uv_cv", fill_lim=daily_rng$uv_cv,
             alpha_lim=NULL,
             xlim=ccz_bbox[1:2], ylim=ccz_bbox[3:4],
             title=paste("", format(date_y, "%b-%d")), out_dim=c(10, 4.5),
             out_f=glue("{dir_figTemp}CCZ_uv_cv_day-{format(date_y, '%m-%d')}.png"),
             tracks=NULL, darkLines=F)
}


daily_df |> 
  filter(!is.na(uv_mn)) |> 
  group_by(lat, lon) |> 
  arrange(desc(uv_mn)) |> 
  slice_head(n=1) |> 
  mutate(yday=yday(date)) |> 
  ggplot() + 
  geom_raster(aes(lon, lat, fill=yday)) + 
  scale_fill_gradientn(colours=cmr$infinity) + 
  geom_sf(data=ccz) + 
  geom_point(data=POI, aes(lon, lat, shape=name), colour="yellow") + 
  theme(legend.position="bottom")




# Images: Eastern focus ---------------------------------------------------

if(reload) {
  cmems_east <- download_cmems(creds=readRDS("data/cmems_cred.rds"),
                               east_bbox, 
                               dates=c(ymd("2023-01-01"), today()),
                               out_nc=glue("data/CMEMS/cmems_east_2023-01-01_present.nc")) |>
    mutate(EKE=pmin(EKE_cm2s2, 3e3))
saveRDS(cmems_east, glue("data/CMEMS/cmems_east_{min(cmems_east$date)}_{max(cmems_east$date)}.rds"))
} else {
  cmems_east <- readRDS(last(dirf("data/CMEMS", "cmems_east.*rds"))) |>
    filter(lat >= 3)
}
if(FALSE) {
  cmems_east |> filter(date==first(date)) |>
    select(lat, lon) |>
    mutate(lat=c(lat), 
           lon=c(lon)) |>
    filter(between(lon, domain_bbox$xmin, domain_bbox$xmax) & 
           between(lat, domain_bbox$ymin, domain_bbox$ymax)) |>
    mutate(id=row_number()) |>
    st_as_sf(coords=c("lon", "lat"), crs=4326) |>
    st_write("data/CMEMS/SMARTEX_domain_centroids.gpkg", append=F)
}



# JC241 mooring longitude
lims <- with(cmems_east |>
               filter(abs(c(lon) + 116.49) < 0.12 &
                      between(lat, 12, 15)) |>
               filter(between(month(date), 2, 7)),
             list(sla=range(sla),
                  uv=range(uv),
                  EKE=range(EKE),
                  u=range(ugos),
                  v=range(vgos),
                  ua=range(ugosa),
                  va=range(vgosa)))

cmems_east |> 
  filter(abs(c(lon) + 116.49) < 0.12, 
         between(lat, 12, 15)) |>
  filter(between(month(date), 2, 7)) |>
  ggplot(aes(date, lat, fill=sla)) + 
  geom_raster() + 
  scale_fill_viridis_c("ssh anomaly (m)") + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_el_2023-FebJul_CMEMS.png", width=10, height=4)

cmems_east |> 
  filter(abs(c(lon) + 116.49) < 0.12, 
         between(lat, 12, 15)) |>
  filter(between(month(date), 2, 7)) |>
  ggplot(aes(date, lat, fill=uv)) + 
  geom_raster() + 
  find_palette("uv", lims=lims$uv) +
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_uv_2023-FebJul_CMEMS.png", width=10, height=4)

cmems_east |> 
  filter(abs(c(lon) + 116.49) < 0.12, 
         between(lat, 12, 15)) |>
  filter(between(month(date), 2, 7)) |>
  ggplot(aes(date, lat, fill=EKE)) + 
  geom_raster() + 
  find_palette("EKE", lims=lims$EKE) +
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_EKE_2023-FebJul_CMEMS.png", width=10, height=4)

cmems_east |> 
  filter(abs(c(lon) + 116.49) < 0.12, 
         between(lat, 12, 15)) |>
  filter(between(month(date), 2, 7)) |>
  ggplot(aes(date, lat, fill=ugos)) + 
  geom_raster() + 
  find_palette("ugos", lims=lims$ugos) +
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_u_2023-FebJul_CMEMS.png", width=10, height=4)

cmems_east |> 
  filter(abs(c(lon) + 116.49) < 0.12, 
         between(lat, 12, 15)) |>
  filter(between(month(date), 2, 7)) |>
  ggplot(aes(date, lat, fill=vgos)) + 
  geom_raster() + 
  find_palette("vgos", lims=lims$vgos) +
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_v_2023-FebJul_CMEMS.png", width=10, height=4)


cmems_east |> 
  filter(abs(c(lon) + 116.49) < 0.12, 
         between(lat, 12, 15)) |>
  filter(between(month(date), 2, 7)) |>
  mutate(uvDir=atan2(vgos, ugos)) |>
  ggplot(aes(date, lat, fill=uvDir)) + 
  geom_raster() + 
  find_palette("uvDir") +
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_uvDir_2023-FebJul_CMEMS.png", width=10, height=4)



# EKE 2D maps
EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
                out_f=c("mean", "sd", "variance"),
                legend=c("EKE mean", "EKE sd", "EKE variance"))
for(k in 1:nrow(EKE_i)) {
  plot_EKE_2D(cmems_east, bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
              east_bbox[1:2], east_bbox[3:4],
              glue("figs/anim/east_EKE_{EKE_i$out_f[k]}_3N-20N.png"), 
              c(10, 5.5))
}
EKE_sf <- cmems_east |>
  group_by(lat, lon) |>
  summarise(EKE_mn=mean(EKE),
            EKE_var=var(EKE),
            EKE_sd=sd(EKE)) |>
  ungroup() |>
  st_as_sf(coords=c("lon", "lat"), crs=4326)
walk(EKE_i$var, ~EKE_sf |> select(all_of(.x), geometry) |>
       stars::st_rasterize() |>
       stars::write_stars(glue("data/CMEMS/{.x}_east_3N-20N.tif")))

timesteps <- sort(unique(cmems_east$date))
cmems_vars <- list("sla"=range(cmems_east$sla, na.rm=T),
                   "EKE"=range(cmems_east$EKE, na.rm=T),
                   # "uv"=range(cmems_east$uv, na.rm=T))
                   "uv"=c(0,1))
tracks_east <- readRDS(last(dirf("data/AVISO", "tracks_east.*rds"))) |>
  filter(type=="anticyclonic")

for(i in seq_along(timesteps)) {
  i_time <- timesteps[i]
  i_df <- cmems_east |> filter(date==i_time)
  i_tracks <- tracks_east |> filter(date <= i_time) |>
    group_by(track) |>
    mutate(active=factor(track==132761 | track==145386, levels=c("FALSE", "TRUE"))) |>
    ungroup()

  for(k in seq_along(cmems_vars)) {
    plot_cmems(df=i_df |> filter(!is.na(sla)), bathy=NULL, land=land, ccz=ccz, POI=POI,
               fill_var=names(cmems_vars)[k], fill_lim=cmems_vars[[k]],
               alpha_lim=NULL,#range(sqrt(abs(cmems_east$sla)), na.rm=T),
               xlim=east_bbox[1:2], ylim=east_bbox[3:4],
               title=i_time, out_dim=c(10, 6),
               out_f=glue("{dir_figTemp}east_{names(cmems_vars)[k]}_{i_time}.png"),
               tracks=i_tracks, darkLines=names(cmems_vars)[k]=="sla")
    gc()
  }
}
rm(cmems_east); gc()




# Images: Full CCZ --------------------------------------------------------

# cmems_ccz <- download_cmems(creds=readRDS("data/cmems_cred.rds"), ccz_bbox) |>
#   mutate(EKE=pmin(EKE_cm2s2, 3e3))
# 
# 
# # EKE 2D maps
# EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
#                 out_f=c("mean", "sd", "variance"),
#                 legend=c("EKE mean", "EKE sd", "EKE variance"))
# for(k in 1:nrow(EKE_i)) {
#   plot_EKE_2D(cmems_ccz, bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
#               ccz_bbox[1:2], ccz_bbox[3:4],
#               glue("figs/anim/FullCCZ_EKE_{EKE_i$out_f[k]}.png"), c(11, 5))
# }
# 
# 
# timesteps <- sort(unique(cmems_ccz$date))
# sla_rng <- range(cmems_ccz$sla, na.rm=T)
# EKE_rng <- range(cmems_ccz$EKE, na.rm=T)
# tracks_ccz <- readRDS(last(dirf("data/AVISO", "tracks_ccz.*rds"))) |>
#   filter(type=="anticyclonic")
# 
# for(i in seq_along(timesteps)) {
#   i_time <- timesteps[i]
#   i_df <- cmems_ccz |> filter(date==i_time)
#   i_tracks <- tracks_ccz |> filter(date <= i_time) |>
#     group_by(track) |>
#     # mutate(active=factor(max(date)==i_time, levels=c("FALSE", "TRUE"))) |>
#     mutate(active=factor(track==132761, levels=c("FALSE", "TRUE"))) |>
#     ungroup()
# 
#   plot_cmems(df=i_df, bathy=bathy, ccz=ccz, POI=POI,
#              fill_var="sla", fill_lim=sla_rng, title=i_time,
#              out_dim=c(10, 4),
#              out_f=glue("{dir_figTemp}CCZ_sla_{i_time}.png"),
#              tracks=i_tracks)
#   gc()
#   plot_cmems(df=i_df, bathy=bathy, ccz=ccz, POI=POI,
#              fill_var="EKE", fill_lim=EKE_rng, title=i_time,
#              out_dim=c(10, 4),
#              out_f=glue("{dir_figTemp}CCZ_EKE_{i_time}.png"),
#              tracks=i_tracks)
#   gc()
# 
# }
# rm(cmems_ccz); gc()



# Images: Eastern focus, multi-year ---------------------------------------

# cmems_east <- download_cmems(creds=readRDS("data/cmems_cred.rds"), east_bbox,
#                              dates=c(ymd("2019-12-01"), today())) |>
#   mutate(EKE=pmin(EKE_cm2s2, 3e3))
# 
# # EKE 2D maps
# EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
#                 out_f=c("mean", "sd", "variance"),
#                 legend=c("EKE mean", "EKE sd", "EKE variance"))
# for(k in 1:nrow(EKE_i)) {
#   plot_EKE_2D(cmems_ccz, bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
#               ccz_bbox[1:2], ccz_bbox[3:4], 
#               glue("figs/anim/FullCCZ_EKE_{EKE_i$out_f[k]}.png"), c(11, 5))
# }
# 
# timesteps <- sort(unique(cmems_east$date))
# sla_rng <- range(cmems_east$sla, na.rm=T)
# EKE_rng <- range(cmems_east$EKE, na.rm=T)
# 
# for(i in seq_along(timesteps)) {
#   i_time <- timesteps[i]
#   i_df <- cmems_east |> filter(date==i_time)
# 
#   plot_cmems(df=i_df, bathy=bathy, ccz=ccz, POI=POI,
#              fill_var="sla", fill_lim=sla_rng,
#              xlim=east_bbox[1:2], ylim=east_bbox[3:4],
#              title=i_time, out_dim=c(10, 4),
#              out_f=glue("{dir_figTemp}eastFull_sla_{i_time}.png"))
#   gc()
#   plot_cmems(df=i_df, bathy=bathy, ccz=ccz, POI=POI,
#              fill_var="EKE", fill_lim=EKE_rng,
#              xlim=east_bbox[1:2], ylim=east_bbox[3:4],
#              title=i_time, out_dim=c(10, 4),
#              out_f=glue("{dir_figTemp}eastFull_EKE_{i_time}.png"))
#   gc()
# 
# }
# rm(cmems_east); gc()




# Images: Full CCZ, multi-year --------------------------------------------

if(reload) {
  cmems_ccz <- download_cmems(creds=readRDS("data/cmems_cred.rds"), ccz_bbox,
                              dates=c(ymd("2019-12-01"), today())) |>
    mutate(EKE=pmin(EKE_cm2s2, 3e3))
  saveRDS(cmems_ccz, glue("data/CMEMS/cmems_ccz_{min(cmems_ccz$date)}_{max(cmems_ccz$date)}.rds"))  
} else {
  cmems_ccz <- readRDS(first(dirf("data/CMEMS", "cmems_ccz_2019-12-01"))) |>
    filter(lat >= 5)
}

# EKE 2D maps
EKE_i <- tibble(var=c("EKE_mn", "EKE_sd", "EKE_var"),
                out_f=c("mean", "sd", "variance"),
                legend=c("EKE mean", "EKE sd", "EKE variance"))
for(k in 1:nrow(EKE_i)) {
  plot_EKE_2D(cmems_ccz, bathy, ccz, POI, EKE_i$var[k], EKE_i$legend[k],
              ccz_bbox[1:2], ccz_bbox[3:4],
              glue("figs/anim/FullCCZ_EKE_{EKE_i$out_f[k]}_5N-23.5N.png"), c(11, 5))
}
EKE_sf <- cmems_ccz |>
  group_by(lat, lon) |>
  summarise(EKE_mn=mean(EKE),
            EKE_var=var(EKE),
            EKE_sd=sd(EKE)) |>
  ungroup() |>
  st_as_sf(coords=c("lon", "lat"), crs=4326)
walk(EKE_i$var, ~EKE_sf |> select(all_of(.x), geometry) |>
       stars::st_rasterize() |>
       stars::write_stars(glue("data/CMEMS/{.x}_CCZ_5N-23.5N.tif")))

# Animations
timesteps <- sort(unique(cmems_ccz$date))
sla_rng <- range(cmems_ccz$sla, na.rm=T)
EKE_rng <- range(cmems_ccz$EKE, na.rm=T)
uv_rng <- c(0,1)

for(i in seq_along(timesteps)) {
  i_time <- timesteps[i]
  i_df <- cmems_ccz |> filter(date==i_time)
  
  plot_cmems(df=i_df |> filter(!is.na(sla)), bathy=bathy, land=land, ccz=ccz, POI=POI,
             fill_var="sla", fill_lim=sla_rng,
             alpha_lim=range(sqrt(abs(cmems_ccz$sla)), na.rm=T),
             title=i_time, out_dim=c(10, 4.75),
             out_f=glue("{dir_figTemp}CCZFull_sla_{i_time}.png"), darkLines=T)
  gc()
  # plot_cmems(df=i_df, bathy=NULL, ccz=ccz, POI=POI,
  #            fill_var="EKE", fill_lim=EKE_rng,  alpha_lim=c(0, max(abs(sla_rng))),
  #            title=i_time,
  #            out_dim=c(10, 4.75),
  #            out_f=glue("{dir_figTemp}CCZFull_EKE_{i_time}.png"), darkLines=F)
  # gc()
  # plot_cmems(df=i_df, bathy=NULL, ccz=ccz, POI=POI,
  #            fill_var="uv", fill_lim=uv_rng,  alpha_lim=c(0, max(abs(sla_rng))),
  #            title=i_time,
  #            out_dim=c(10, 4.75),
  #            out_f=glue("{dir_figTemp}CCZFull_uv_{i_time}.png"), darkLines=F)
  # gc()

}
rm(cmems_ccz); gc()



# Animate -----------------------------------------------------------------

img_path <- dir_figTemp
out_path <- "figs/anim/"

fps <- 20
sets <- c(outer(c("CCZ_", "east_", "CCZFull_", "eastFull_"), 
                c("sla", "EKE", "uv"), 
                paste0))
# sets <- c(outer(c("CCZ_EKE_mean_", 
#                   "CCZ_uv_mean_",
#                   "CCZ_EKE_sdInterannual_", 
#                   "CCZ_uv_sdInterannual_",
#                   "CCZ_uv_cv_",
#                   "CCZ_EKE_cv_"), 
#                 c("year", "month", "week", "day"), 
#                 paste0))
for(i in sets) {
  # fps <- switch(str_split_fixed(i, "_", 4)[1,4],
  #               "year"=3,
  #               "month"=3,
  #               "week"=8,
  #               "day"=50)
  tic()
  try({
    av_encode_video(dirf(img_path, i), 
                    glue("figs/anim/{i}.mp4"), 
                    framerate=fps)
  })
  toc()
  gc()
}
