# HYCOM download
# SMARTEX
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(ncdf4)
library(sf)
library(stars)
theme_set(theme_bw() + theme(panel.grid=element_blank()))

dirf("code/fn/") |> walk(source)

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
date_seq <- seq(ymd("2020-01-01"), ymd("2021-12-31"), by=1) 

thredds_url <- "http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0"
hy_nc <- nc_open(thredds_url)
hy_i <- list(lon=ncvar_get(hy_nc, "lon"),
             lat=ncvar_get(hy_nc, "lat"),
             time=ncvar_get(hy_nc, "time"),
             depth=ncvar_get(hy_nc, "depth"))
nc_close(hy_nc)




# download daily nc -------------------------------------------------------

for(i in seq_along(date_seq)) {
  gc()
  d <- date_seq[i]
  ymdh_i <- as.numeric(difftime(as_datetime(paste(d, "00:00:00")),
                                as_datetime("2000-01-01 00:00:00"),
                                units="hours"))
  cat("Starting", as.character(d), "at", format(Sys.time()), "\n")

  file_base <- "hycom_glby_930_SMARTEX"

  download_hycom(
    url_base=thredds_url,
    hy_i=hy_i,
    bbox=domain_bbox,
    time_range=ymdh_i + c(0, 21),
    depth_range=c(1, 40),
    out_nc=glue("{dirs$nc}{file_base}_{d}.nc")
  )
}




# process to R dataframes -------------------------------------------------

for(i in seq_along(date_seq)) {
  gc()
  d <- date_seq[i]
  ymdh_i <- as.numeric(difftime(as_datetime(paste(d, "00:00:00")), 
                                as_datetime("2000-01-01 00:00:00"),
                                units="hours"))
  file_base <- "hycom_glby_930_SMARTEX"
  
  download_hycom(
    hy_i=hy_i,
    bbox=domain_bbox,
    time_range=ymdh_i + c(0, 21), 
    depth_range=c(1, 40),
    out_nc=glue("{dirs$nc}{file_base}_{d}.nc"), 
    out_rds=glue("{dirs$rds}{file_base}_{d}.rds"), 
    out_rng=glue("{dirs$rds}{file_base}_ranges_{d}.rds")
  )
}




# extract JC241 mooring ---------------------------------------------------

jc241 <- dirf(dirs$rds, "930_SMARTEX.*3D.rds") |>
  grep("ranges", x=_, value=T, invert=T) |>
  map_dfr(~readRDS(.x) |>
            filter(id==2404) |> # minimized abs(latDiff), abs(lonDiff)
            select(time, depth, lon, lat, u, v, temperature, salinity)) |>
  filter(depth != 5000) |>
  group_by(depth) |>
  mutate(v=if_else(v > 28, mean(v), v)) |> # 2020-06-20 09:00 -- ??
  ungroup()
write_csv(jc241, "hycom_JC241_mooring_2020-23.csv")

jc241 |>
  mutate(depth=paste0("d_", depth)) |>
  select(time, depth, u) |>
  filter(year(time)==2023 & between(month(time), 3, 6)) |>
  pivot_wider(names_from=depth, values_from=u) |>
  select(-time) |>
  write_csv("jc241_u_focal.csv")
jc241 |>
  mutate(depth=paste0("d_", depth)) |>
  select(time, depth, v) |>
  filter(year(time)==2023 & between(month(time), 3, 6)) |>
  pivot_wider(names_from=depth, values_from=v) |>
  select(-time) |>
  write_csv("jc241_v_focal.csv")


bgr <- dirf(dirs$rds, "930_SMARTEX.*3D.rds") |>
  grep("ranges", x=_, value=T, invert=T) |>
  map_dfr(~readRDS(.x) |>
            filter(id==1901) |>
            select(time, depth, lon, lat, u, v, temperature, salinity))
write_csv(bgr, "hycom_bgr_mooring__2020-23.csv")












# TEMP: Plots, etc --------------------------------------------------------

jc241 <- read_csv("hycom_JC241_mooring_2020-23.csv") |>  # 2020-06-20 09:00 -- ??
  group_by(depth) |>
  mutate(v=if_else(v > 28, mean(v), v)) |>
  ungroup()

eddy_aviso_dates <- tibble(xmin=ymd("2023-03-28"), 
                           xmax=ymd("2023-04-26"),
                           ymin=-Inf, ymax=Inf)

mag_scale <- 4
p <- jc241 |>
  mutate(depth=as.numeric(factor(depth)),
         uvDir=atan2(v,u)) |>
  ggplot() +
  geom_rect(data=eddy_aviso_dates, 
            aes(xmin=as_datetime(xmin), xmax=as_datetime(xmax), ymin=ymin, ymax=ymax), 
            alpha=0.1, colour=NA) +
  geom_segment(aes(x=time, xend=time+u*mag_scale*3600*24, 
                   y=-depth, yend=-depth+v*mag_scale, colour=uvDir),
               arrow=arrow(length=unit(0.02, "inches")), linewidth=0.5) +
  scale_colour_gradientn(
    name="Current\ndirection",
    colours=readRDS("data/cmr_cmaps.RDS")$infinity, limits=c(-pi, pi),
    breaks=c(-pi, -pi/2, 0, pi/2, pi),
    labels=c("W", "S", "E", "N", "W")) +
  scale_y_continuous("Depth (m)", breaks=-(1:39), labels=unique(jc241$depth)) +
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.2),
        legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5))
ggsave("figs/JC241_arrows_3-hourly_2020-23.png", p, width=60, height=9, limitsize=F)

mag_scale <- 6
jc241 |>
  mutate(date=date(time)) |>
  group_by(date, depth) |>
  summarise(u=mean(u), v=mean(v)) |>
  ungroup() |>
  mutate(depth=as.numeric(factor(depth)),
         uvDir=atan2(v,u)) |>
  ggplot() +
  geom_rect(data=eddy_aviso_dates, 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            alpha=0.1, colour=NA) +
  geom_segment(aes(x=date, xend=date+u*mag_scale, 
                   y=-depth, yend=-depth+v*mag_scale, colour=uvDir),
               arrow=arrow(length=unit(0.05, "inches")), linewidth=0.75) +
  scale_colour_gradientn(
    name="Current\ndirection",
    colours=readRDS("data/cmr_cmaps.RDS")$infinity, limits=c(-pi, pi),
    breaks=c(-pi, -pi/2, 0, pi/2, pi),
    labels=c("W", "S", "E", "N", "W")) +
  scale_y_continuous("Depth (m)", breaks=-(1:39), labels=unique(jc241$depth)) +
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.2),
        legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5))
ggsave("figs/JC241_arrows_daily_2020-23.png", p, width=80, height=9, limitsize=F)


jc241 |>
  mutate(depth=as.numeric(factor(depth)),
         uvDir=atan2(v,u),
         uv=sqrt(u^2 + v^2)) |>
  ggplot(aes(time, -depth, fill=uvDir)) + geom_raster() +
  scale_fill_gradientn(
    name="Current\ndirection",
    colours=readRDS("data/cmr_cmaps.RDS")$infinity, limits=c(-pi, pi),
    breaks=c(-pi, -pi/2, 0, pi/2, pi),
    labels=c("W", "S", "E", "N", "W")) +
  scale_y_continuous("Depth (m)", breaks=-(1:39), labels=unique(jc241$depth)) +
  scale_x_datetime(date_breaks="14 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.2),
        legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5))
ggsave("figs/JC241_uvDir_2020-23.png", width=29, height=9)


jc241 |>
  mutate(uv=sqrt(u^2 + v^2)) |>
  select(time, depth, u, v, uv) |>
  pivot_longer(3:5, names_to="Component", values_to="Velocity") |>
  mutate(Component=factor(Component, levels=c("u", "v", "uv"))) |>
  ggplot() + 
  geom_rect(data=eddy_aviso_dates, 
            aes(xmin=as_datetime(xmin), xmax=as_datetime(xmax), ymin=ymin, ymax=ymax), 
            alpha=0.1, colour=NA) +
  geom_hline(yintercept=0) +
  geom_line(aes(time, Velocity, group=depth, colour=log1p(depth))) + 
  scale_colour_viridis_c("Depth (m)", direction=-1, 
                         breaks=log1p(c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)),
                         labels=c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)) + 
  scale_x_datetime(date_breaks="14 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  labs(y="m/s") +
  facet_grid(Component~., scales="free_y") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241_velocities-3h_2020-23.png", width=20, height=10)

jc241 |>
  mutate(uv=sqrt(u^2 + v^2)) |>
  select(time, depth, u, v, uv) |>
  pivot_longer(3:5, names_to="Component", values_to="Velocity") |>
  mutate(Component=factor(Component, levels=c("u", "v", "uv"))) |>
  group_by(depth, Component) |>
  mutate(Velocity=c(scale(Velocity))) |>
  ungroup() |>
  ggplot() + 
  geom_rect(data=eddy_aviso_dates, 
            aes(xmin=as_datetime(xmin), xmax=as_datetime(xmax), ymin=ymin, ymax=ymax), 
            alpha=0.1, colour=NA) +
  geom_hline(yintercept=0) +
  geom_line(aes(time, Velocity, group=depth, colour=log1p(depth))) + 
  scale_colour_viridis_c("Depth (m)", direction=-1, 
                         breaks=log1p(c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)),
                         labels=c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)) + 
  scale_x_datetime(date_breaks="14 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  labs(y="z(Velocity) by depth") +
  facet_grid(Component~., scales="free_y") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241_z-velocities-3h_2020-23.png", width=20, height=10)



jc241 |>
  select(time, depth, u, v) |>
  mutate(date=date(time)) |>
  group_by(date, depth) |> 
  summarise(u=mean(u), v=mean(v), uv=sqrt(u^2 + v^2)) |>
  ungroup() |>
  pivot_longer(3:5, names_to="Component", values_to="Velocity") |>
  mutate(Component=factor(Component, levels=c("u", "v", "uv"))) |>
  ggplot() + 
  geom_rect(data=eddy_aviso_dates, 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            alpha=0.1, colour=NA) +
  geom_hline(yintercept=0) +
  geom_line(aes(date, Velocity, group=depth, colour=log1p(depth))) + 
  scale_colour_viridis_c("Depth (m)", direction=-1, 
                         breaks=log1p(c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)),
                         labels=c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)) + 
  scale_x_date(date_breaks="14 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  labs(y="Daily mean (m/s)") +
  facet_grid(Component~., scales="free_y") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241_velocities-daily_2020-23.png", width=20, height=10)


jc241 |>
  select(time, depth, u, v) |>
  mutate(date=date(time)) |>
  group_by(date, depth) |> 
  summarise(u=mean(u), v=mean(v), uv=sqrt(u^2 + v^2)) |>
  ungroup() |>
  pivot_longer(3:5, names_to="Component", values_to="Velocity") |>
  mutate(Component=factor(Component, levels=c("u", "v", "uv"))) |>
  group_by(depth, Component) |>
  mutate(Velocity=c(scale(Velocity))) |>
  ungroup() |>
  ggplot() + 
  geom_rect(data=eddy_aviso_dates, 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            alpha=0.1, colour=NA) +
  geom_hline(yintercept=0) +
  geom_line(aes(date, Velocity, group=depth, colour=log1p(depth))) + 
  scale_colour_viridis_c("Depth (m)", direction=-1, 
                         breaks=log1p(c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)),
                         labels=c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)) + 
  scale_x_date(date_breaks="14 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  labs(y="z(daily mean velocity) by depth") +
  facet_grid(Component~., scales="free_y") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241_z-velocities-daily_2020-23.png", width=20, height=10)

jc241 |>
  select(time, depth, u, v) |>
  group_by(depth) |>
  mutate(u=zoo::rollmean(u, k=7*8, fill=NA),
         v=zoo::rollmean(v, k=7*8, fill=NA),
         uv=sqrt(u^2 + v^2)) |>
  pivot_longer(3:5, names_to="Component", values_to="Velocity") |>
  mutate(Component=factor(Component, levels=c("u", "v", "uv"))) |>
  ggplot() + 
  geom_rect(data=eddy_aviso_dates, 
            aes(xmin=as_datetime(xmin), xmax=as_datetime(xmax), ymin=ymin, ymax=ymax), 
            alpha=0.1, colour=NA) +
  geom_hline(yintercept=0) +
  geom_line(aes(time, Velocity, group=depth, colour=log1p(depth))) + 
  scale_colour_viridis_c("Depth (m)", direction=-1, 
                         breaks=log1p(c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)),
                         labels=c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)) + 
  scale_x_datetime(date_breaks="14 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  labs(y="Rolling weekly mean velocity by depth") +
  facet_grid(Component~., scales="free_y") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241_velocities-weekly_2020-23.png", width=20, height=10)

jc241 |>
  select(time, depth, u, v) |>
  group_by(depth) |>
  mutate(u=zoo::rollmean(u, k=7*8, fill=NA),
         v=zoo::rollmean(v, k=7*8, fill=NA),
         uv=sqrt(u^2 + v^2)) |>
  pivot_longer(3:5, names_to="Component", values_to="Velocity") |>
  mutate(Component=factor(Component, levels=c("u", "v", "uv"))) |>
  group_by(depth, Component) |>
  mutate(Velocity=c(scale(Velocity))) |>
  ungroup() |>
  ggplot() + 
  geom_rect(data=eddy_aviso_dates, 
            aes(xmin=as_datetime(xmin), xmax=as_datetime(xmax), ymin=ymin, ymax=ymax), 
            alpha=0.1, colour=NA) +
  geom_hline(yintercept=0) +
  geom_line(aes(time, Velocity, group=depth, colour=log1p(depth))) + 
  scale_colour_viridis_c("Depth (m)", direction=-1, 
                         breaks=log1p(c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)),
                         labels=c(0, 10, 25, 50, 100, 250, 500, 1000, 4000)) + 
  scale_x_datetime(date_breaks="14 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  labs(y="z(rolling weekly mean velocity) by depth") +
  facet_grid(Component~., scales="free_y") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241_z-velocities-weekly_2020-23.png", width=20, height=10)





depth_i <- tibble(depth=unique(jc241$depth)) |>
  mutate(depth0=depth - (depth-lag(depth, default=0))/2,
         depth1=depth - (depth-lead(depth, default=0))/2)

jc241 |>
  select(time, depth, u, v) |>
  pivot_longer(3:4, names_to="Component", values_to="Velocity") |>
  mutate(Component=factor(Component, levels=c("u", "v"))) |>
  left_join(depth_i) |>
  ggplot(aes(xmin=time-1.5*60*60, xmax=time+1.5*60*60, ymin=depth0, ymax=depth1, fill=Velocity)) + 
  geom_rect(colour=NA) + 
  scale_fill_gradient2("m/s") + 
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  scale_y_reverse() +
  facet_grid(Component~., scales="free_y") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))

jc241 |>
  mutate(uv=sqrt(u^2 + v^2)) |>
  select(time, depth, u, v, uv) |>
  pivot_longer(3:5, names_to="Component", values_to="Velocity") |>
  mutate(Component=factor(Component, levels=c("u", "v", "uv"))) |>
  group_by(depth, Component) |>
  mutate(Velocity=c(scale(Velocity))) |>
  ungroup() |>
  left_join(depth_i) |>
  ggplot(aes(xmin=time-1.5*60*60, xmax=time+1.5*60*60, ymin=depth0, ymax=depth1, fill=Velocity)) + 
  geom_rect(colour=NA) + 
  scale_fill_viridis_c("sd", direction=-1) + 
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  scale_y_reverse() +
  facet_grid(Component~., scales="free_y") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))






# corrplot ----------------------------------------------------------------

# u 
jc241 |> 
  filter(depth != 5000) |> 
  select(time, depth, u) |> 
  pivot_wider(names_from=depth, values_from=u) |> 
  select(-time) |> 
  cor() |> 
  corrplot::corrplot(title="u, 3-hourly")

jc241 |> 
  mutate(date=as_date(time)) |> 
  filter(depth != 5000) |> 
  group_by(date, depth) |>
  summarise(u=mean(u)) |> 
  ungroup() |>
  select(date, depth, u) |> 
  pivot_wider(names_from=depth, values_from=u) |> 
  select(-date) |> 
  cor() |> 
  corrplot::corrplot(title="u, daily means")

# v 
jc241 |> 
  filter(depth != 5000) |> 
  select(time, depth, v) |> 
  pivot_wider(names_from=depth, values_from=v) |> 
  select(-time) |> 
  cor() |> 
  corrplot::corrplot(title="v, 3-hourly")

jc241 |> 
  mutate(date=as_date(time)) |> 
  filter(depth != 5000) |> 
  group_by(date, depth) |>
  summarise(v=mean(v)) |> 
  ungroup() |>
  select(date, depth, v) |> 
  pivot_wider(names_from=depth, values_from=v) |> 
  select(-date) |> 
  cor() |> 
  corrplot::corrplot(title="v, daily means")

# uv 
jc241 |> 
  filter(depth != 5000) |> 
  mutate(uv=sqrt(u^2+v^2)) |>
  select(time, depth, uv) |> 
  pivot_wider(names_from=depth, values_from=uv) |> 
  select(-time) |> 
  cor() |> 
  corrplot::corrplot(title="uv, 3-hourly")

jc241 |> 
  mutate(date=as_date(time)) |> 
  filter(depth != 5000) |> 
  group_by(date, depth) |>
  summarise(uv=mean(sqrt(u^2+v^2))) |> 
  ungroup() |>
  select(date, depth, uv) |> 
  pivot_wider(names_from=depth, values_from=uv) |> 
  select(-date) |> 
  cor() |> 
  corrplot::corrplot(title="uv, daily means")


# # Better to do in matlab with mscohere() and cpsd()
# 
# d0 <- jc241 |> filter(depth==0)
# d4000 <- jc241 |> filter(depth==4000)
# # Convert velocities to complex numbers
# d0_complex <- complex(real=d0$u, imaginary=d0$v)
# d4000_complex <- complex(real=d4000$u, imaginary=d4000$v)
# 
# library(freqdom)
# spec_dens <- spectral.density(d0_complex, d4000_complex)
# plot(spec_dens$freq, spec_dens$operators[1,1,])
# 
# spec_dens2 <- spectrum(cbind(d0_complex, d4000_complex))







# # Perform Fourier transform
# d0_complex_spectrum <- fft(d0_complex)
# d4000_complex_spectrum <- fft(d4000_complex)
# 
# # Calculate cross-spectrum
# d0_cross_spectrum <- Conj(d0_complex_spectrum) * d0_complex_spectrum
# d4000_cross_spectrum <- Conj(d4000_complex_spectrum) * d4000_complex_spectrum
# d0_d4000_cross <- Conj(d0_complex_spectrum) * d4000_complex_spectrum
# 
# # Calculate coherence magnitude
# coherence_magnitude <- abs(d0_d4000_cross) / (abs(d0_cross_spectrum) * abs(d4000_cross_spectrum))
# 
# # Calculate coherence phase
# coherence_phase <- Arg(d0_d4000_cross)
# 
# plot(coherence_magnitude, type="b")
# plot(coherence_phase, type="b")








# surface patterns --------------------------------------------------------

surf_df <- full_join(
  dirf(dirs$rds, "SMARTEX.*2D.rds") |>
    grep("ranges", x=_, value=T, invert=T) |>
    grep("2023-0[2-7]", x=_, value=T) |>
    map_dfr(~readRDS(.x) |>
              select(time, lat, lon, id, surf_el)),
  dirf(dirs$rds, "SMARTEX.*3D.rds") |>
    grep("ranges", x=_, value=T, invert=T) |>
    grep("2023-0[2-7]", x=_, value=T) |>
    map_dfr(~readRDS(.x) |>
              filter(depth==0) |>
              select(time, lat, lon, id, u, v))
) |>
  mutate(uv=sqrt(u^2 + v^2))

limits <- list(
  u=range(filter(surf_df, abs(lon-243.52)<0.01)$u),
  v=range(filter(surf_df, abs(lon-243.52)<0.01)$v),
  uv=range(filter(surf_df, abs(lon-243.52)<0.01)$uv)
)

library(colorspace)
surf_df %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(time, lat, fill=u)) + 
  geom_raster() + 
  find_palette("u", limits$u) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_u_2023-FebJul_HYCOM.png", width=10, height=4)

surf_df %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(time, lat, fill=v)) + 
  geom_raster() + 
  find_palette("v", limits$v) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_v_2023-FebJul_HYCOM.png", width=10, height=4)

surf_df %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(time, lat, fill=atan2(v,u))) + 
  geom_raster() + 
  find_palette("uvDir") + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_uvDir_2023-FebJul_HYCOM.png", width=10, height=4)

surf_df %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(time, lat, fill=sqrt(u^2+v^2))) + 
  geom_raster() + 
  scale_fill_viridis_c(option="turbo", begin=0.1, end=0.9) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_uv_2023-FebJul_HYCOM.png", width=10, height=4)

surf_df %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(time, lat, fill=surf_el)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_el_2023-FebJul_HYCOM.png", width=10, height=4)




hy_daily <- surf_df |>
  mutate(date=date(time)) |>
  group_by(date, lat, lon, id) |>
  summarise(across(everything(), mean)) |>
  ungroup()
limits <- list(
  u=range(filter(hy_daily, abs(lon-243.52)<0.01)$u),
  v=range(filter(hy_daily, abs(lon-243.52)<0.01)$v),
  uv=range(filter(hy_daily, abs(lon-243.52)<0.01)$uv)
)

library(colorspace)
hy_daily %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(date, lat, fill=u)) + 
  geom_raster() + 
  find_palette("u", limits$u) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_u_2023-FebJul_HYCOM-daily.png", width=10, height=4)

hy_daily %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(date, lat, fill=v)) + 
  geom_raster() + 
  find_palette("v", limits$v) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_v_2023-FebJul_HYCOM-daily.png", width=10, height=4)

hy_daily %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(date, lat, fill=atan2(v,u))) + 
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
ggsave("figs/JC241-lon_surf_uvDir_2023-FebJul_HYCOM-daily.png", width=10, height=4)

hy_daily %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(date, lat, fill=sqrt(u^2+v^2))) + 
  geom_raster() + 
  scale_fill_viridis_c(option="turbo", begin=0.1, end=0.9) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_uv_2023-FebJul_HYCOM-daily.png", width=10, height=4)

hy_daily %>% 
  filter(abs(lon-243.52) < 0.01) |> 
  ggplot(aes(date, lat, fill=surf_el)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_el_2023-FebJul_HYCOM-daily.png", width=10, height=4)



daily_latSlice <- bind_rows(
  surf_df %>% 
    filter(abs(lon-243.52) < 0.01) |>
    mutate(date=date(time)) |>
    group_by(date, lat) |>
    summarise(u=mean(u), 
              v=mean(v)) |>
    ungroup() |>
    mutate(source="HYCOM"),
  readRDS(last(dirf("data/CMEMS", "cmems_east.*rds"))) |>
    filter(abs(c(lon) + 116.49) < 0.12 &
           between(lat, 11.7, 15.2)) |>
    filter(between(month(date), 2, 7)) |> 
    select(date, lat, ugos, vgos) |>
    rename(u=ugos, v=vgos) |>
    mutate(source="CMEMS") 
) |>
  mutate(uv=sqrt(u^2 + v^2),
         uvDir=atan2(v, u))

limits <- with(daily_latSlice, 
               list(u=range(u), 
                    v=range(v),
                    uv=range(uv)))

ggplot(daily_latSlice, aes(date, lat, fill=u)) + 
  geom_raster() + 
  find_palette("u", limits$u) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  facet_grid(source~.) + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_u_2023-FebJul_daily.png", width=10, height=8)

ggplot(daily_latSlice, aes(date, lat, fill=v)) + 
  geom_raster() + 
  find_palette("v", limits$v) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  facet_grid(source~.) + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_v_2023-FebJul_daily.png", width=10, height=8)

ggplot(daily_latSlice, aes(date, lat, fill=uv)) + 
  geom_raster() + 
  find_palette("uv", limits$uv) + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  facet_grid(source~.) + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_uv_2023-FebJul_daily.png", width=10, height=8)

ggplot(daily_latSlice, aes(date, lat, fill=uvDir)) + 
  geom_raster() + 
  find_palette("uvDir") + 
  geom_hline(yintercept=13.88, linetype=3) + 
  labs(x="Date", y="Latitude (deg. N)", title="JC241 mooring longitude") + 
  scale_x_date(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  facet_grid(source~.) + 
  theme(legend.position="bottom", 
        legend.key.height=unit(0.25, "cm"), 
        legend.key.width=unit(3, "cm"),
        axis.text.x=element_text(angle=270, hjust=0.5, vjust=0.5),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
ggsave("figs/JC241-lon_surf_uvDir_2023-FebJul_daily.png", width=10, height=8)





library(colorspace)
figs_temp <- "D:/SMARTEX/figs/temp/"
ccz <- st_read("data/CCZ_areas/ccz_outline.gpkg")

depth_subsets <- c(0)
hycom_par <- tibble(par=c("u", "v", "uv", "uvDir", "surf_el"),
                    flab=c("u", "v", "uv", "uvDir", "surf_el"),
                    dark=c(T, T, T, T, T))

surf_df <- surf_df |>
  mutate(EKE_std=uv,
         uvDir=atan2(v, u), 
         depth=0,
         lon=lon-360) 
timesteps <- unique(surf_df$time)
rng <- with(surf_df, 
            list("u"=range(u, na.rm=T),
                 "v"=range(v, na.rm=T),
                 "uv"=range(uv, na.rm=T),
                 "surf_el"=range(surf_el, na.rm=T)))
lon_diff <- mean(diff(sort(unique(surf_df$lon))))
timesteps <- sort(unique(surf_df$time))
lon_rng <- range(surf_df$lon)
lat_rng <- range(surf_df$lat)
for(i in seq_along(timesteps)) {
  gc()
  i_date <- as_date(timesteps)[i]
  i_hour <- hour(timesteps)[i]
  i_df <- filter(surf_df, date(time)==i_date & hour(time)==i_hour)
  
  for(k in 1:nrow(hycom_par)) {
    plot_hycom_depthPanels(
      df=i_df, depths=depth_subsets, ccz=ccz, POI=POI,
      fill_var=hycom_par$par[k], fill_lim=rng[[hycom_par$par[k]]],
      xlim=lon_rng, ylim=lat_rng,
      title=paste0(i_date, " ", str_pad(i_hour, 2, "left", "0"), ":00"), 
      out_dim=c(8, 4.5),
      out_f=glue("{figs_temp}hycomSurface_{hycom_par$flab[k]}_",
                 "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
      tracks=NULL, darkLines=T)
  }
  plot_hycom_depthPanels(
    df=i_df, depths=depth_subsets, ccz=ccz, POI=POI,
    fill_var="uvDir", alpha_lim=rng$uv,
    xlim=lon_rng, ylim=lat_rng,
    title=paste0(i_date, " ", str_pad(i_hour, 2, "left", "0"), ":00"), 
    out_dim=c(8, 4.5),
    out_f=glue("{figs_temp}hycomSurface_uvDir_",
               "{i_date}_{str_pad(i_hour, 2, 'left', '0')}.png"),
    tracks=NULL, darkLines=T)
}


library(tictoc); library(av); library(sevcheck); library(glue)
img_path <- "figs/anim/temp"
img_path <- "D:/SMARTEX/figs/temp/"
out_path <- "figs/anim_temp/"

fps <- 20
sets <- c(outer(c("hycomSurface_"),#, "eddy4000m_", "eddy2000m_", "eddy1000m_", "eddy500m_"), 
                c("u", "v", "u2v2", "uv", "surf_el",
                  "z-u", "z-v", "z-u2v2", "z-uv",
                  "uvDir", "EKE", "z-EKE"), 
                paste0))

for(i in sets) {
  try({
    cat("starting", i, "\n")
    tic()
    av::av_encode_video(dirf(img_path, paste0(i, "_")), 
                        glue("{out_path}SMARTEX_{i}.mp4"),
                        framerate=fps)
    toc()
  }, silent=T)
  gc()
}


