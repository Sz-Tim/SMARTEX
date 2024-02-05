# ARGO floats
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


# jc257 bbox: 9 to 18N, -118 to -106E


# setup
library(tidyverse); library(glue); library(ncdf4); library(sf); library(gsw)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
theme_set(theme_bw() + theme(panel.grid=element_blank()))
walk(dirf("code/fn", ".R$"), source)


domain_bbox <- list(xmin=-119, xmax=-114, ymin=12, ymax=15)
domain_rect <- cbind(unlist(domain_bbox[1:2])[c(1,2,2,1,1)], 
                     unlist(domain_bbox[3:4])[c(1,1,2,2,1)]) |>
  list() |> st_polygon() |> st_sfc(crs=4326) |> st_sf()
st_geometry(domain_rect) <- "geometry"
POI <- read_csv("data/SMARTEX_locations_copy.csv") |>
  filter(grepl("BGR|Long|Atoll|001|002|006", Name_Full))
ccz_east <- st_read("data/CCZ_areas/ccz_outline.gpkg") |>
  st_crop(c(xmin=-125, xmax=-90, ymin=7, ymax=20))
eez_east <- st_read("data/CCZ_areas/CCZ_EEZs.shp") |>
  st_crop(c(xmin=-125, xmax=-90, ymin=7, ymax=20))
jc257_eddy <- st_read("data/cruise2024/cruise2024_eddies_SSH_peak.gpkg") |>
  mutate(dataset="jc257 (Dec-)")
jc257_eddy_now <- st_read("data/cruise2024/cruise2024_eddies_radius.gpkg") |> 
  slice_tail(n=1) |>
  mutate(dataset="jc257 (Dec-)")



# load data ---------------------------------------------------------------

argo_f <- dirrf("data/ARGO/aoml", "*.nc$")
# argo_f <- dirrf("data/ARGO/", "*.nc$")

argo_df <- map_dfr(argo_f, process_argo_nc) |>
  mutate(date=ymd(mtime),
         absSal=gsw_SA_from_SP(if_else(is.na(S), Si, S), if_else(is.na(P), Pi, P), lon, lat),
         consTemp=gsw_CT_from_pt(absSal, if_else(is.na(T), Ti, T)),
         rho=gsw_rho(absSal, consTemp, if_else(is.na(P), Pi, P)),
         sigma0=gsw_sigma0(absSal, consTemp),
         depth=-gsw_z_from_p(if_else(is.na(P), Pi, P), lat))

traj_df <- argo_df |>
  select(-c("Pi", "Ti", "Si", "P", "T", "S", "absSal", "consTemp", "rho", "sigma0", "depth", "obs")) |>
  group_by(PN, FSN, dir_id, lat, lon, mtime) |>
  slice_head(n=1) |>
  ungroup() |>
  arrange(dir_id, mtime) |>
  mutate(date=ymd(mtime),
         dataset=if_else(mtime < ymd("2023-09-01"), 
                         "2023 eddy (Feb-Aug)", 
                         "jc257 (Dec-)"),
         month=month(mtime)) |>
  group_by(dir_id) |>
  mutate(latest=row_number() == n()) |>
  ungroup() |>
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=F) %>%
  mutate(outside_domain=!map_dbl(st_within(., domain_rect, sparse=T), length),
         km_from_domain=as.numeric(c(st_distance(., domain_rect)))/1e3 * outside_domain) |>
  st_drop_geometry()

# 2023 Eddy
argo_2023eddy_i <- traj_df |>
  filter(dataset=="2023 eddy (Feb-Aug)") |>
  filter(between(lon, domain_bbox$xmin-3, domain_bbox$xmax+3),
         between(lat, domain_bbox$ymin-3, domain_bbox$ymax+3)) |>
  mutate(fnm=str_remove(fnm, "data/ARGO/")) |>
  select(dir_id, fnm, PN, CN, FSN, juld, mtime, lon, lat, profile, km_from_domain) 
argo_2023eddy_i |>
  write_csv("data/ARGO/eddy2023/ARGO_metadata_2023-eddy.csv")

# copy files to new directory
if(FALSE) {
  for(i in 1:nrow(argo_2023eddy_i)) {
    dir_i <- dirname(glue("data/ARGO/{argo_2023eddy_i$fnm[i]}"))
    dir_new <- str_replace(dir_i, "ARGO/aoml", "ARGO/eddy2023/aoml")
    dir.create(dir_new, recursive=T, showWarnings=F)
    file.copy(glue("data/ARGO/{argo_2023eddy_i$fnm[i]}"), dir_new)
  }
}

# jc257
jc257_i <- traj_df |>
  filter(dataset=="jc257 (Dec-)") |>
  filter(between(lon, -118, -105), between(lat, 10, 16)) |>
  mutate(fnm=str_remove(fnm, "data/ARGO/")) |>
  select(dir_id, fnm, PN, CN, FSN, juld, mtime, lon, lat, profile) 
jc257_i |>
  write_csv("data/ARGO/jc257/ARGO_metadata_jc257.csv")

# copy files to new directory
if(FALSE) {
  for(i in 1:nrow(jc257_i)) {
    dir_i <- dirname(glue("data/ARGO/{jc257_i$fnm[i]}"))
    dir_new <- str_replace(dir_i, "ARGO/aoml", "ARGO/jc257/aoml")
    dir.create(dir_new, recursive=T, showWarnings=F)
    file.copy(glue("data/ARGO/{jc257_i$fnm[i]}"), dir_new)
  }
}





# plots -------------------------------------------------------------------

traj_df |>
  ggplot() + 
  geom_sf(data=domain_rect, alpha=0.25, colour=NA, fill="grey30") +
  geom_sf(data=ccz_east) +
  geom_sf(data=eez_east, linewidth=1, fill=NA) +
  geom_sf(data=jc257_eddy, colour="red", shape=1, size=0.5) +
  geom_sf(data=jc257_eddy_now, colour="red", fill=NA) +
  geom_point(aes(lon, lat, colour=month)) + 
  geom_path(aes(lon, lat, group=dir_id)) + 
  geom_text(data=traj_df |> filter(latest), 
            aes(lon, lat, label=dir_id), size=3, hjust=0, nudge_x=0.2, colour="red") +
  facet_wrap(~dataset) +
  scale_colour_viridis_c(option="turbo", breaks=1:12) +
  xlim(-125, -104) + ylim(9, 18) +
  theme(legend.position="bottom",
        legend.key.width=unit(2, 'cm'),
        legend.key.height=unit(0.2, 'cm'))
ggsave("figs/ARGO_map.png", width=9, height=4)



traj_df |>
  filter(dataset=="jc257 (Dec-)") |>
  filter(between(lon, -118, -105), between(lat, 10, 16)) |>
  ggplot() + 
  geom_sf(data=ccz_east, colour="grey") +
  geom_sf(data=eez_east, linewidth=1, fill=NA, colour="grey") +
  geom_sf(data=jc257_eddy, colour="steelblue", shape=1, size=0.5) +
  geom_sf(data=jc257_eddy_now, colour="steelblue", fill=NA) +
  geom_point(aes(lon, lat, colour=latest)) + 
  geom_path(aes(lon, lat, group=dir_id)) + 
  geom_text(aes(lon, lat, label=dir_id, size=latest),
            hjust=0, nudge_x=0.2, colour="red3") +
  facet_wrap(~dataset) +
  scale_colour_manual(values=c("FALSE"="grey30", "TRUE"="red3"), guide="none") +
  scale_size_manual(values=c("FALSE"=0.001, "TRUE"=2.5), guide="none") +
  xlim(-115, -104) + ylim(10, 16) +
  theme(legend.position="bottom",
        legend.key.width=unit(2, 'cm'),
        legend.key.height=unit(0.2, 'cm'),
        axis.title=element_blank(),
        panel.grid=element_line(colour="grey85", linewidth=0.1))
ggsave("data/ARGO/jc257/ARGO_jc257.png", width=7, height=4)




traj_df |>
  filter(dataset=="2023 eddy (Feb-Aug)") |>
  filter(between(lon, domain_bbox$xmin-3, domain_bbox$xmax+3),
         between(lat, domain_bbox$ymin-3, domain_bbox$ymax+3)) |>
  ggplot() + 
  geom_sf(data=domain_rect, alpha=0.25, colour=NA, fill="grey30") +
  geom_sf(data=ccz_east) +
  geom_sf(data=eez_east, linewidth=1, fill=NA) +
  geom_point(aes(lon, lat, colour=km_from_domain)) + 
  geom_path(aes(lon, lat, group=dir_id)) + 
  scale_colour_viridis_c(option="turbo") +
  xlim(domain_bbox$xmin-3, domain_bbox$xmax+3) + ylim(domain_bbox$ymin-3, domain_bbox$ymax+3) +
  theme(legend.position="bottom",
        legend.key.width=unit(1.5, 'cm'),
        legend.key.height=unit(0.2, 'cm'),
        axis.title=element_blank(),
        panel.grid=element_line(colour="grey85", linewidth=0.1))
ggsave("figs/ARGO_2023.png", width=5, height=5)



traj_df |> 
  filter(between(lon, -120, -115),
         between(lat, 12, 15)) |>
  ggplot() + 
  geom_rect(data=as_tibble(domain_bbox), alpha=0.25,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
  geom_sf(data=ccz_east) +
  geom_point(aes(lon, lat, colour=dir_id)) + 
  geom_path(aes(lon, lat, group=dir_id)) +
  theme(legend.position=c(0.8, 0.8))
ggsave("figs/ARGO_2023_eddy_mooring.png", width=4, height=4)


argo_2023_eddy <- argo_df |>
  filter(dir_id %in% c("R5906797"),
         between(mtime, ymd("2023-03-01"), ymd("2023-06-30"))) |>
  mutate(yday=yday(mtime),
         month=month(mtime, label=TRUE),
         date=factor(mtime))

argo_2023_eddy_i <- argo_2023_eddy |>
  group_by(mtime) |>
  slice_head(n=1) |>
  ungroup()

hycom_df <- map_dfr(1:nrow(argo_2023_eddy_i), 
                    ~dirf("data/HYCOM/", glue("SMARTEX_{argo_2023_eddy_i$mtime[.x]}_3D.rds")) |>
                      readRDS() |>
                      mutate(lon=lon-360) |>
                      mutate(lon_diff=abs(lon - argo_2023_eddy_i$lon[.x]),
                             lat_diff=abs(lat - argo_2023_eddy_i$lat[.x])) |>
                      filter(lon_diff==min(lon_diff),
                             lat_diff==min(lat_diff))) |>
  mutate(mtime=as_date(time)) |>
  mutate(pressure=gsw_p_from_z(-depth, lat),
         absSal=gsw_SA_from_SP(salinity, pressure, lon, lat),
         consTemp=gsw_CT_from_pt(absSal, temperature),
         rho=gsw_rho(absSal, consTemp, pressure),
         sigma0=gsw_sigma0(absSal, consTemp)) |>
  group_by(mtime, depth) |>
  summarise(across(where(is.numeric), mean)) |>
  ungroup() |>
  mutate(yday=yday(mtime),
         month=month(mtime, label=TRUE),
         date=factor(mtime)) 

depth_lim <- 250
p_ls <- vector("list", 6)
p_ls[[1]] <-  argo_2023_eddy |>
  filter(P <= depth_lim) |>
  ggplot() + 
  geom_point(aes(consTemp, P, colour=date), size=0.5) + 
  geom_path(aes(consTemp, P, group=paste(mtime, profile), colour=date)) + 
  geom_point(data=hycom_df |> filter(pressure <= depth_lim), 
             aes(consTemp, pressure, colour=date), shape=1, alpha=0.5) + 
  scale_colour_viridis_d(option="turbo") +
  scale_y_reverse() +
  guides(colour=guide_legend(nrow=2)) +
  labs(x="Conservative temperature (C)", y="Pressure (db)") +
  theme(panel.grid.major=element_line(colour="grey90"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.key.width=unit(1, 'cm'),
        legend.key.height=unit(0.2, 'cm'))
p_ls[[2]] <- argo_2023_eddy |>
  filter(P <= depth_lim) |>
  ggplot() + 
  geom_point(aes(absSal, P, colour=date), size=0.5) + 
  geom_path(aes(absSal, P, group=paste(mtime, profile), colour=date)) + 
  geom_point(data=hycom_df |> filter(pressure <= depth_lim), 
             aes(absSal, pressure, colour=date), shape=1, alpha=0.5) + 
  scale_colour_viridis_d(option="turbo") +
  scale_y_reverse() +
  guides(colour=guide_legend(nrow=2)) +
  labs(x="Absolute salinity (g/kg)", y="Pressure (db)") +
  theme(panel.grid.major=element_line(colour="grey90"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.key.width=unit(1, 'cm'),
        legend.key.height=unit(0.2, 'cm'))
p_ls[[3]] <- argo_2023_eddy |>
  filter(P <= depth_lim) |>
  ggplot() + 
  geom_point(aes(absSal, consTemp, colour=date), size=0.5) + 
  geom_path(aes(absSal, consTemp, group=paste(mtime, profile), colour=date)) + 
  geom_point(data=hycom_df |> filter(pressure <= depth_lim), 
             aes(absSal, consTemp, colour=date), shape=1, alpha=0.5) + 
  scale_colour_viridis_d(option="turbo") +
  guides(colour=guide_legend(nrow=2)) +
  labs(x="Absolute salinity (g/kg)", y="Conservative temperature (C)") +
  theme(panel.grid.major=element_line(colour="grey90"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.key.width=unit(1, 'cm'),
        legend.key.height=unit(0.2, 'cm'))
p_ls[[4]] <- argo_2023_eddy |>
  filter(P <= depth_lim) |>
  ggplot() + 
  geom_point(aes(rho, P, colour=date), size=0.5) + 
  geom_path(aes(rho, P, group=paste(mtime, profile), colour=date)) + 
  geom_point(data=hycom_df |> filter(pressure <= depth_lim), 
             aes(rho, pressure, colour=date), shape=1, alpha=0.5) + 
  scale_colour_viridis_d(option="turbo") +
  scale_y_reverse() +
  guides(colour=guide_legend(nrow=2)) +
  labs(x="in-situ density rho (kg/m3)", y="Pressure (db)") +
  theme(panel.grid.major=element_line(colour="grey90"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.key.width=unit(1, 'cm'),
        legend.key.height=unit(0.2, 'cm'))
p_ls[[5]] <- argo_2023_eddy |>
  filter(P <= depth_lim) |>
  ggplot() + 
  geom_point(aes(sigma0, P, colour=date), size=0.5) + 
  geom_path(aes(sigma0, P, group=paste(mtime, profile), colour=date)) + 
  geom_point(data=hycom_df |> filter(pressure <= depth_lim), 
             aes(sigma0, pressure, colour=date), shape=1, alpha=0.5) + 
  scale_colour_viridis_d(option="turbo") +
  scale_y_reverse() +
  guides(colour=guide_legend(nrow=2)) +
  labs(x="sigma0 (kg/m3)", y="Pressure (db)") +
  theme(panel.grid.major=element_line(colour="grey90"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.key.width=unit(1, 'cm'),
        legend.key.height=unit(0.2, 'cm'))
p_ls[[6]] <- argo_2023_eddy_i |>
  ggplot() + 
  geom_point(data=POI |> filter(grepl("Long mooring", Name_Full)),
             aes(lon, lat), shape=4, size=3) +
  geom_sf(data=ccz_east |> st_crop(c(xmin=-118.25, xmax=-114, ymin=11.5, ymax=16))) +
  geom_path(aes(lon, lat)) + 
  geom_point(aes(lon, lat, colour=date)) + 
  geom_point(data=hycom_df |> group_by(date) |> slice_head(n=1) |> ungroup(),
             aes(lon, lat, colour=date), size=3, shape=1) + 
  scale_colour_viridis_d(option="turbo") +
  guides(colour=guide_legend(nrow=2, title=element_blank())) +
  theme(axis.title=element_blank(), 
        legend.title=element_blank(),
        legend.position="bottom",
        legend.key.width=unit(1, 'cm'),
        legend.key.height=unit(0.2, 'cm'))
ggpubr::ggarrange(plotlist=p_ls, ncol=3, nrow=2, 
                  common.legend=T, legend="bottom")
ggsave(glue("figs/ARGO_2023_mooring_{first(argo_2023_eddy$dir_id)}.png"), width=9, height=6)






# cmems -------------------------------------------------------------------

if(FALSE) {
  library(ncdf4)
  east_bbox <- list(xmin=-130, xmax=-90, ymin=2, ymax=20)
  cmems_2023 <- download_cmems(creds=readRDS("data/cmems_cred.rds"),
                               east_bbox, 
                               dates=c(ymd("2023-01-01"), 
                                       ymd("2023-12-31")),
                               out_nc=glue("data/CMEMS/cmems_east_2023.nc"))
  saveRDS(cmems_2023, "data/CMEMS/cmems_east_2023.rds")
}
cmems_2023 <- readRDS("data/CMEMS/cmems_east_2023.rds") |>
  filter(between(lon, domain_bbox$xmin-6, domain_bbox$xmax+6),
         between(lat, domain_bbox$ymin-5, domain_bbox$ymax+6)) 
subset_2023 <- argo_2023eddy_i |> mutate(date=ymd(mtime)) |> filter(km_from_domain < 200)
dates_2023 <- sort(unique(subset_2023$date))

depth_lim <- 200
argo_df_sub <- argo_df |>
  mutate(fnm=str_remove(fnm, "data/ARGO/")) |>
  filter(fnm %in% subset_2023$fnm)
traj_df_sub <- traj_df |>
  mutate(fnm=str_remove(fnm, "data/ARGO/")) |>
  filter(fnm %in% subset_2023$fnm)
rng <- list(sla=range(cmems_2023$sla, na.rm=T),
            uv=range(cmems_2023$uv, na.rm=T),
            absSal=range(filter(argo_df_sub, P <= depth_lim, )$absSal, na.rm=T),
            consTemp=range(filter(argo_df_sub, P <= depth_lim)$consTemp, na.rm=T),
            rho=range(filter(argo_df_sub, P <= depth_lim)$rho, na.rm=T),
            sigma0=range(filter(argo_df_sub, P <= depth_lim)$sigma0, na.rm=T),
            P=rev(pmin(range(filter(argo_df_sub, P <= depth_lim)$P, na.rm=T), depth_lim)))

for(i in seq_along(dates_2023)) {
  argo_i <- argo_df_sub |> filter(ymd(mtime)==dates_2023[i])
  traj_i <- traj_df_sub |> filter(ymd(mtime)==dates_2023[i]) 
  cmems_i <- cmems_2023 |> filter(date==ymd(dates_2023[i]))
  
  if(nrow(argo_i) > 0) {
    p <- ggplot(traj_i) +
      geom_raster(data=cmems_i, aes(lon, lat, fill=sla)) + 
      find_palette("sla", lims=rng$sla) + 
      geom_sf(data=ccz_east, colour="grey40", linewidth=0.2) +
      geom_sf(data=domain_rect, colour="grey40", fill=NA, linewidth=0.4) + 
      geom_point(aes(lon, lat), size=2) + 
      geom_text(aes(lon, lat, label=dir_id), hjust=0, nudge_x=0.25) +
      xlim(domain_bbox$xmin-6, domain_bbox$xmax+6) + 
      ylim(domain_bbox$ymin-5, domain_bbox$ymax+5) +
      ggtitle(dates_2023[i]) +
      theme(legend.position="bottom",
            legend.key.width=unit(1.5, 'cm'),
            legend.key.height=unit(0.2, 'cm'),
            axis.title=element_blank(),
            panel.grid=element_line(colour="grey85", linewidth=0.1))
    ggsave(glue("figs/ARGO/sla_{dates_2023[i]}.png"), p, width=5, height=5)
    
    p <- ggplot(traj_i) +
      geom_raster(data=cmems_i, aes(lon, lat, fill=uv)) + 
      find_palette("uv", lims=rng$uv) + 
      geom_sf(data=ccz_east, colour="grey90", linewidth=0.2) +
      geom_sf(data=domain_rect, colour="grey90", fill=NA, linewidth=0.4) + 
      geom_point(aes(lon, lat), size=2) + 
      geom_text(aes(lon, lat, label=dir_id), hjust=0, nudge_x=0.25) +
      xlim(domain_bbox$xmin-6, domain_bbox$xmax+6) + 
      ylim(domain_bbox$ymin-5, domain_bbox$ymax+5) +
      ggtitle(dates_2023[i]) +
      theme(legend.position="bottom",
            legend.key.width=unit(1.5, 'cm'),
            legend.key.height=unit(0.2, 'cm'),
            axis.title=element_blank(),
            panel.grid=element_line(colour="grey85", linewidth=0.1))
    ggsave(glue("figs/ARGO/uv_{dates_2023[i]}.png"), p, width=5, height=5)
    
    
    # profiles
    p_ls <- vector("list", 6)
    p_ls[[1]] <-  argo_i |>
      filter(P <= depth_lim) |>
      ggplot() + 
      geom_point(aes(consTemp, P, colour=dir_id), size=0.5) + 
      geom_path(aes(consTemp, P, group=paste(fnm, profile), colour=dir_id)) + 
      scale_colour_viridis_d(option="turbo") +
      scale_y_reverse(limits=rng$P) + xlim(rng$consTemp) +
      guides(colour=guide_legend(nrow=2)) +
      labs(x="Conservative temperature (C)", y="Pressure (db)") +
      theme(panel.grid.major=element_line(colour="grey90"),
            legend.title=element_blank(),
            legend.position="bottom",
            legend.key.width=unit(1, 'cm'),
            legend.key.height=unit(0.2, 'cm'))
    p_ls[[2]] <- argo_i |>
      filter(P <= depth_lim) |>
      ggplot() + 
      geom_point(aes(absSal, P, colour=dir_id), size=0.5) + 
      geom_path(aes(absSal, P, group=paste(fnm, profile), colour=dir_id)) + 
      scale_colour_viridis_d(option="turbo") +
      scale_y_reverse(limits=rng$P) + xlim(rng$absSal) +
      guides(colour=guide_legend(nrow=2)) +
      labs(x="Absolute salinity (g/kg)", y="Pressure (db)") +
      theme(panel.grid.major=element_line(colour="grey90"),
            legend.title=element_blank(),
            legend.position="bottom",
            legend.key.width=unit(1, 'cm'),
            legend.key.height=unit(0.2, 'cm'))
    p_ls[[3]] <- argo_i |>
      filter(P <= depth_lim) |>
      ggplot() + 
      geom_point(aes(absSal, consTemp, colour=dir_id), size=0.5) + 
      geom_path(aes(absSal, consTemp, group=paste(fnm, profile), colour=dir_id)) + 
      scale_colour_viridis_d(option="turbo") +
      guides(colour=guide_legend(nrow=2)) +
      ylim(rng$consTemp) + xlim(rng$absSal) +
      labs(x="Absolute salinity (g/kg)", y="Conservative temperature (C)") +
      theme(panel.grid.major=element_line(colour="grey90"),
            legend.title=element_blank(),
            legend.position="bottom",
            legend.key.width=unit(1, 'cm'),
            legend.key.height=unit(0.2, 'cm'))
    p_ls[[4]] <- argo_i |>
      filter(P <= depth_lim) |>
      ggplot() + 
      geom_point(aes(rho, P, colour=dir_id), size=0.5) + 
      geom_path(aes(rho, P, group=paste(fnm, profile), colour=dir_id)) + 
      scale_colour_viridis_d(option="turbo") +
      scale_y_reverse(limits=rng$P) + xlim(rng$rho) +
      guides(colour=guide_legend(nrow=2)) +
      labs(x="in-situ density rho (kg/m3)", y="Pressure (db)") +
      theme(panel.grid.major=element_line(colour="grey90"),
            legend.title=element_blank(),
            legend.position="bottom",
            legend.key.width=unit(1, 'cm'),
            legend.key.height=unit(0.2, 'cm'))
    p_ls[[5]] <- argo_i |>
      filter(P <= depth_lim) |>
      ggplot() + 
      geom_point(aes(sigma0, P, colour=dir_id), size=0.5) + 
      geom_path(aes(sigma0, P, group=paste(fnm, profile), colour=dir_id)) + 
      scale_colour_viridis_d(option="turbo") +
      scale_y_reverse(limits=rng$P) + xlim(rng$sigma0) +
      guides(colour=guide_legend(nrow=2)) +
      labs(x="sigma0 (kg/m3)", y="Pressure (db)") +
      theme(panel.grid.major=element_line(colour="grey90"),
            legend.title=element_blank(),
            legend.position="bottom",
            legend.key.width=unit(1, 'cm'),
            legend.key.height=unit(0.2, 'cm'))
    p_ls[[6]] <- traj_i |>
      ggplot() + 
      geom_raster(data=cmems_i, aes(lon, lat, fill=sla)) + 
      find_palette("sla", lims=rng$sla) + 
      geom_sf(data=ccz_east, colour="grey30", linewidth=0.2) +
      geom_sf(data=domain_rect, colour="grey30", fill=NA, linewidth=0.4) + 
      xlim(domain_bbox$xmin-6, domain_bbox$xmax+6) + 
      ylim(domain_bbox$ymin-5, domain_bbox$ymax+5) +
      ggtitle(dates_2023[i]) +
      geom_point(aes(lon, lat, colour=dir_id)) + 
      geom_text(aes(lon, lat, label=dir_id), hjust=0, nudge_x=0.25, size=2) +
      scale_colour_viridis_d(option="turbo") +
      guides(colour=guide_legend(nrow=2, title=element_blank())) +
      theme(axis.title=element_blank(), 
            legend.title=element_blank(),
            legend.position="bottom",
            legend.key.width=unit(1, 'cm'),
            legend.key.height=unit(0.2, 'cm'))
    p <- ggpubr::ggarrange(plotlist=p_ls, ncol=3, nrow=2, 
                           common.legend=T, legend="bottom")
    ggsave(glue("figs/ARGO/ARGO_profiles_{dates_2023[i]}.png"), p, width=9, height=6)
  }
  gc()
}
