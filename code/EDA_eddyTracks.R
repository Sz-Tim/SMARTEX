# Eddy tracks
# SMARTEX
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


# Eddy tracks from AVISO+

# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(terra)
library(sf)
library(ncdf4)
library(ggnewscale)
library(colorspace)
library(magick)
library(tictoc)

theme_set(theme_bw() + theme(panel.grid=element_blank()))

dirf("code/fn/") |> walk(source)


POI <- read_csv("data/SMARTEX_locations_copy.csv") |>
  filter(grepl("BGR|Long|Atoll|001|002", Name_Full))

ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")
east_bbox <- list(xmin=-130, xmax=-90, ymin=2, ymax=20)
ccz_bbox <- list(xmin=-160, xmax=-90, ymin=0, ymax=23.5)

bathy_reclass <- cbind(c(-Inf, seq(-6000, 0, by=500)),
                       c(seq(-6000, 0, by=500), Inf),
                       seq(-6.5, 0, by=0.5)*1e3)
bathy <- rast("data/bathymetry/gebco_2023_n23.5_s0.0_w-160.0_e-90.0.tif") |>
  aggregate(fact=20) |>
  classify(bathy_reclass) |>
  as.polygons() |> 
  st_as_sf() |> 
  rename_with(~"elev", .cols=1)



# CCZ: 2022-10-01 to today()
ccz_tracks <- get_eddyTracks(ccz_bbox, dates=c(ymd("2022-10-01"), today())) |>
  group_by(track) |>
  mutate(nObs=diff(range(observation_number)),
         minLon=min(lon),
         maxLon=max(lon)) |>
  ungroup() |>
  filter(nObs > (16*7)) 
saveRDS(ccz_tracks,
        glue("data/AVISO/tracks_ccz_{min(ccz_tracks$date)}_{max(ccz_tracks$date)}.rds"))

ccz_tracks |> filter(type=="anticyclonic") |>
  ggplot() + 
  # geom_sf(data=bathy, aes(fill=elev), colour=NA, alpha=0.5) +
  # scale_fill_distiller(type="seq", palette="Greys", direction=1, guide="none") +
  # new_scale_fill() +
  geom_path(aes(lon, lat, colour=amplitude, group=track), linewidth=1, alpha=0.8) + 
  geom_sf(data=ccz) + 
  geom_point(data=POI, aes(lon, lat, shape=name), size=4) +
  # scale_shape_manual("", values=c(1,4)) + 
  scale_colour_viridis_c(option="plasma") +
  xlim(ccz_bbox[[1]], ccz_bbox[[2]]) + ylim(ccz_bbox[[3]], ccz_bbox[[4]])


# East focus: 2022-10-01 to today()
east_tracks <- get_eddyTracks(east_bbox, dates=c(ymd("2022-10-01"), today())) |>
  group_by(track) |>
  mutate(nObs=diff(range(observation_number)),
         minLon=min(lon),
         maxLon=max(lon)) |>
  ungroup() |>
  filter(nObs > (16*7)) 
saveRDS(east_tracks,
        glue("data/AVISO/tracks_east_{min(east_tracks$date)}_{max(east_tracks$date)}.rds"))
east_tracks |> filter(type=="anticyclonic") |>
  # filter(track==132761) |>
  ggplot() + 
  # geom_sf(data=bathy, aes(fill=elev), colour=NA, alpha=0.5) +
  # scale_fill_distiller(type="seq", palette="Greys", direction=1, guide="none") +
  # new_scale_fill() +
  # geom_point(aes(lon, lat, colour=amplitude, size=effective_radius/1e3), shape=1) +
  # geom_linerange(aes(lon, ymin=lat-effective_radius/111e3, ymax=lat+effective_radius/111e3,
                     # colour=amplitude), alpha=0.2, linewidth=1) +
  # geom_ribbon(aes(lon, ymin=lat-effective_radius/111e3, ymax=lat+effective_radius/111e3,
  #                    group=track), alpha=0.3, colour=NA, fill=scales::muted("blue")) +
  geom_path(aes(lon, lat, colour=amplitude, group=track), linewidth=1, alpha=0.8) + 
  geom_sf(data=ccz) + 
  geom_point(data=POI, aes(lon, lat, shape=name), size=4) +
  # scale_shape_manual("", values=c(1,4)) + 
  scale_colour_viridis_c(option="plasma") +
  scale_size_area("Effective area (km2)", max_size=10) +
  xlim(east_bbox[[1]], east_bbox[[2]]) + ylim(east_bbox[[3]], east_bbox[[4]])


domain_bbox <- list(xmin=-119, xmax=-114, ymin=12, ymax=15)
domain_tracks <- get_eddyTracks(domain_bbox, dates=c(ymd("2022-10-01"), today())) |>
  group_by(track) |>
  mutate(nObs=diff(range(observation_number)),
         minLon=min(lon),
         maxLon=max(lon)) |>
  ungroup()
  # filter(nObs > (16*7)) 
saveRDS(domain_tracks,
        glue("data/AVISO/tracks_POI_{min(domain_tracks$date)}_{max(domain_tracks$date)}.rds"))

# track 145386

tracks_sf <- domain_tracks |>
  mutate(lon=c(lon), lat=c(lat), effective_radius=c(effective_radius)) |>
  st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326) |>
  st_buffer(dist=domain_tracks$effective_radius)

tracks_sf |> 
  filter(year(date)==2023) |>
  ggplot() + 
  geom_sf(aes(colour=type), fill=NA) +
  geom_sf(data=ccz) + 
  geom_point(data=POI, aes(lon, lat, shape=name), size=4) +
  facet_wrap(~month(date, label=T), nrow=2) + 
  xlim(domain_bbox[[1]], domain_bbox[[2]]) + ylim(domain_bbox[[3]]-2, domain_bbox[[4]]+1)
ggsave("figs/AVISO_eddies_2023.png", width=14, height=6)


tracks_sf |> 
  filter(month(date)==7 & year(date)==2023) |>
  ggplot() + 
  geom_sf(aes(colour=day(date), linetype=type), fill=NA, linewidth=0.6) +
  scale_colour_viridis_c(option="turbo") +
  scale_linetype_manual(values=c(1,3)) +
  geom_sf(data=ccz) + 
  geom_point(data=POI, aes(lon, lat, shape=name), size=4) +
  ggtitle("AVISO eddies: July 2023") +
  xlim(domain_bbox[[1]], domain_bbox[[2]]) + ylim(domain_bbox[[3]]-2, domain_bbox[[4]]+1)
ggsave("figs/AVISO_eddies_2023-July.png", width=6, height=6)

tracks_sf |> 
  filter(month(date)==4 & year(date)==2023) |>
  ggplot() + 
  geom_sf(aes(colour=day(date), linetype=type), fill=NA, linewidth=0.6) +
  scale_colour_viridis_c(option="turbo") +
  scale_linetype_manual(values=c(1,3)) +
  geom_sf(data=ccz) + 
  geom_point(data=POI, aes(lon, lat, shape=name), size=4) +
  ggtitle("AVISO eddies: April 2023") +
  xlim(domain_bbox[[1]], domain_bbox[[2]]) + ylim(domain_bbox[[3]]-2, domain_bbox[[4]]+1)
ggsave("figs/AVISO_eddies_2023-April.png", width=6, height=6)

