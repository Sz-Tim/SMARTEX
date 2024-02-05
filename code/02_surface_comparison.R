# Surface: HYCOM vs. CMEMS
# SMARTEX
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


# setup
library(tidyverse); library(glue); library(sf)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
theme_set(theme_bw() + theme(panel.grid=element_blank()))

dirf("code/fn/") |> walk(source)

POI <- read_csv("data/SMARTEX_locations_copy.csv") |>
  filter(grepl("BGR|Long|Atoll|001|002|006", Name_Full))
figs_temp <- "D:/SMARTEX/figs/temp/"
ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")

domain_bbox <- list(xmin=-119, xmax=-114, ymin=12, ymax=15)

# dirs <- list(
#   nc="D:/SMARTEX/data/HYCOM/",
#   rds="D:/SMARTEX/data/HYCOM/",
#   rng="D:/SMARTEX/data/HYCOM/"
# )
dirs <- list(
  nc="data/HYCOM/",
  rds="data/HYCOM/",
  rng="data/HYCOM/"
)


bottom_hycom <- full_join(
  dirf(dirs$rds, "SMARTEX_2023-0[1-3].*_2D.rds") |>
    grep("ranges", x=_, value=T, invert=T) |>
    map_dfr(~readRDS(.x) |>
              select(time, lat, lon, id, ends_with("bottom"))),
  dirf(dirs$rds, "SMARTEX_2023-0[1-3].*_3D.rds") |>
    grep("ranges", x=_, value=T, invert=T) |>
    map_dfr(~readRDS(.x) |>
              filter(!is.na(u)) |>
              arrange(depth) |>
              group_by(id) |>
              slice_tail(n=1) |>
              ungroup() |>
              select(time, lat, lon, id, depth, u, v))
) |>
  mutate(uv=sqrt(u^2 + v^2))


surf_hycom <- full_join(
  dirf(dirs$rds, "SMARTEX_.*_2D.rds") |>
    grep("ranges", x=_, value=T, invert=T) |>
    map_dfr(~readRDS(.x) |>
              select(time, lat, lon, id, surf_el)),
  dirf(dirs$rds, "SMARTEX_.*_3D.rds") |>
    grep("ranges", x=_, value=T, invert=T) |>
    map_dfr(~readRDS(.x) |>
              filter(depth==0) |>
              select(time, lat, lon, id, u, v))
) |>
  mutate(uv=sqrt(u^2 + v^2))
gc()
saveRDS(surf_hycom, "D:/SMARTEX/data/surface_hycom.rds")
rm(surf_hycom)
gc()

if(FALSE) {
  library(ncdf4)
  cmems_SMARTEX <- download_cmems(creds=readRDS("data/cmems_cred.rds"),
                                  map2(domain_bbox, list(-1, 1, -1, 1), ~.x + .y), 
                               dates=c(ymd("2020-01-01"), 
                                       ymd("2023-12-31")),
                               out_nc=glue("data/CMEMS/cmems_SMARTEX_2020-2023.nc"))
  saveRDS(cmems_SMARTEX, "data/CMEMS/CMEMS_SMARTEX.rds")
  SMARTEX_CMEMS_centroids <- cmems_SMARTEX |>
    filter(date==first(date)) |>
    select(lat, lon) |>
    mutate(lat=c(lat), lon=c(lon)) |>
    mutate(CMEMS_id=row_number()) |>
    st_as_sf(coords=c("lon", "lat"), crs=4326)
  st_write(SMARTEX_CMEMS_centroids, "data/CMEMS/SMARTEX_domain_CMEMS_centroids.gpkg")
}
if(FALSE) {
  cmems_jc241 <- cmems_SMARTEX |>
    mutate(lonDiff=abs(lon - POI$lon[2]), 
           latDiff=abs(lat - POI$lat[2])) |>
    arrange(latDiff, lonDiff) |>
    group_by(date) |>
    slice_head(n=1) |>
    ungroup()
  cmems_jc241 |>
    select(ugos, vgos) |>
    write_csv("jc241_cmems_2020-2023.csv")
  cmems_jc241 |>
    filter(year(date)==2023 & between(month(date), 3, 6)) |>
    select(ugos, vgos) |>
    write_csv("jc241_cmems_focal.csv")
}

  

SMARTEX_CMEMS_grid <- st_read("data/CMEMS/SMARTEX_domain_CMEMS_grid.gpkg") %>%
  mutate(lon=round(st_coordinates(st_centroid(.))[,1], 3),
         lat=round(st_coordinates(st_centroid(.))[,2], 3))
surf_CMEMS <- dirf("data/CMEMS", "SMARTEX.*rds") |> last() |> readRDS() |>
  select(date, lat, lon, sla, ugos, vgos, uv) |>
  # filter(year(date)==2023 & between(month(date), 2, 7)) |>
  filter(!is.na(uv)) |>
  mutate(lat=c(lat), lon=c(lon)) |>
  inner_join(SMARTEX_CMEMS_grid |> st_drop_geometry(), by=c("lon", "lat"))

surf_hycom <- readRDS("D:/SMARTEX/data/surface_hycom.rds") |>
  mutate(time=c(time), 
         lat=c(lat),
         lon=c(lon)-360)
hycom_CMEMS_i <- surf_hycom |>
  filter(time==first(time)) |>
  select(lon, lat, id) |>
  st_as_sf(coords=c("lon", "lat"), crs=4326) |>
  st_join(SMARTEX_CMEMS_grid |> select(-lon, -lat), join=st_nearest_feature) |>
  st_drop_geometry()
  



# Summarize HYCOM to CMEMS ------------------------------------------------

hycom_daily_agg <- surf_hycom |>
  right_join(hycom_CMEMS_i, by="id") |>
  mutate(date=date(time)) |>
  group_by(CMEMS_id, date) |>
  summarise(u=mean(u),
            v=mean(v),
            surf_el=mean(surf_el)) |>
  ungroup() |>
  mutate(uv=sqrt(u^2 + v^2),
         source="HYCOM") |>
  group_by(CMEMS_id) |>
  mutate(sla=surf_el - mean(surf_el)) |>
  ungroup() |>
  select(-surf_el) |>
  bind_rows(surf_CMEMS |> 
              select(CMEMS_id, date, ugos, vgos, uv, sla) |>
              rename(u=ugos, v=vgos) |>
              mutate(source="CMEMS"))

hycom_daily_agg |>
  mutate(uv=sqrt(u^2 + v^2)) |>
  pivot_longer(c("u", "v", "uv", "sla"), names_to="variable", values_to="value") |>
  pivot_wider(names_from="source", values_from="value") |>
  left_join(SMARTEX_CMEMS_grid |> st_drop_geometry()) |>
  select(lon, lat, date, variable, CMEMS, HYCOM) |>
  write_csv("SMARTEX_surface_hycom-vs-cmems.csv")



# Scatterplots ------------------------------------------------------------

hycom_daily_agg |>
  pivot_longer(c("u", "v", "uv", "sla"), names_to="variable", values_to="value") |>
  pivot_wider(names_from="source", values_from="value") |>
  mutate(month=month(date, label=T)) |>
  left_join(SMARTEX_CMEMS_grid |> st_drop_geometry()) |>
  ggplot(aes(HYCOM, CMEMS)) + 
  geom_point(alpha=0.1, shape=1) + 
  geom_abline() +
  facet_grid(variable~month) + 
  xlim(-1, 1) + ylim(-1, 1) +
  coord_equal()

hycom_daily_agg |>
  pivot_longer(c("u", "v", "uv", "sla"), names_to="variable", values_to="value") |>
  pivot_wider(names_from="source", values_from="value") |>
  mutate(month=month(date, label=T)) |>
  left_join(SMARTEX_CMEMS_grid |> st_drop_geometry()) |>
  ggplot(aes(HYCOM, CMEMS)) + 
  geom_density2d(contour_var="ndensity") + 
  geom_abline() +
  facet_grid(variable~month) + 
  xlim(-1, 1) + ylim(-1, 1) +
  coord_equal()

hycom_daily_agg |>
  pivot_longer(c("u", "v", "uv", "sla"), names_to="variable", values_to="value") |>
  pivot_wider(names_from="source", values_from="value") |>
  left_join(SMARTEX_CMEMS_grid |> st_drop_geometry()) |>
  mutate(month=month(date, label=T)) |>
  ggplot(aes(HYCOM, CMEMS, colour=month)) + 
  geom_point(shape=1, alpha=0.5) + 
  geom_abline() +
  scale_colour_viridis_d(option="turbo", begin=0.1, end=0.9) +
  guides(colour=guide_legend(override.aes=list(shape=19, alpha=1, size=1.5))) +
  facet_wrap(~variable) + 
  xlim(-1, 1) + ylim(-1, 1) +
  coord_equal()




# Correlation -------------------------------------------------------------

hycom_daily_agg |>
  mutate(uv=sqrt(u^2 + v^2)) |>
  pivot_longer(c("u", "v", "uv", "sla"), names_to="variable", values_to="value") |>
  pivot_wider(names_from="source", values_from="value") |>
  group_by(CMEMS_id, variable) |>
  summarise(r=cor(HYCOM, CMEMS, use="pairwise")) |>
  ungroup() |>
  left_join(SMARTEX_CMEMS_grid |> st_drop_geometry()) |>
  ggplot(aes(lon, lat, fill=r)) + 
  geom_raster() + 
  geom_point(data=POI |> filter(grepl("mooring", name)), fill=NA) + 
  scale_fill_gradient2() +
  facet_wrap(~variable)





# Animations --------------------------------------------------------------

hycom_cmems_diff <- hycom_daily_agg |>
  group_by(CMEMS_id, date) |>
  summarise(across(all_of(c("u", "v", "uv", "sla")), 
                   ~first(.x)-last(.x), 
                   .names="{.col}_diff")) |>
  ungroup() |>
  left_join(SMARTEX_CMEMS_grid |> st_drop_geometry()) |>
  rename(el_diff=sla_diff)
  
lims <- with(hycom_cmems_diff, 
             list(u_diff=range(u_diff),
                  v_diff=range(v_diff),
                  uv_diff=range(uv_diff),
                  el_diff=range(el_diff)))

date_seq <- unique(hycom_cmems_diff$date)
for(i in seq_along(date_seq)) {
  d <- date_seq[i]
  i_df <- filter(hycom_cmems_diff, date==d)
  for(k in seq_along(lims)) {
    plot_cmems(i_df, bathy=NULL, land=NULL, ccz, POI, 
               fill_var=names(lims)[k], fill_lim=lims[[k]], 
               title=glue("HYCOM - CMEMS: {d}  {names(lims)[k]}"), 
               xlim=range(hycom_cmems_diff$lon), ylim=range(hycom_cmems_diff$lat), 
               out_f=glue("{figs_temp}{names(lims)[k]}_{d}.png"), 
               out_dim=c(10, 5.5),
               tracks=NULL, darkLines=TRUE)
  }
  
}


library(tictoc); library(av); library(sevcheck); library(glue)
img_path <- "figs/anim/temp"
img_path <- "D:/SMARTEX/figs/temp/"
out_path <- "figs/anim_temp/"

fps <- 5
sets <- c("u_diff", "v_diff", "uv_diff", "el_diff")

for(i in sets) {
  try({
    cat("starting", i, "\n")
    tic()
    av::av_encode_video(dirf(img_path, paste0("^", i, "_")), 
                        glue("{out_path}{i}.mp4"),
                        framerate=fps)
    toc()
  }, silent=T)
  gc()
}
