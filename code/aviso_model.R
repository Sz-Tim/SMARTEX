## -----------------------------------------------------------------------------
## project
## scriptDescription
## Tim Szewczyk
## -----------------------------------------------------------------------------

library(tidyverse); library(brms); library(sf); library(sevcheck); library(ncdf4); library(glue)
dirf("code/fn/") |> walk(source)
east_bbox <- list(xmin=-130, xmax=-90, ymin=7, ymax=20)

d <- 1
tracks_df <- get_eddyTracks(east_bbox, dates=c(ymd("2020-10-01"), today())) |> 
  filter(month(date) %in% c(10:12, 1:5)) |>
  filter(type=="anticyclonic") |>
  filter(effective_radius > 100000) |>
  group_by(track) |>
  mutate(N=n()) |>
  ungroup() |>
  filter(N > 60) |> 
  select(track, lon, lat, date, effective_radius) |>
  mutate(track=paste0("aviso_", track)) |>
  mutate(across(everything(), c)) |>
  rename(rad=effective_radius) |>
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |>
  st_transform(8859) %>% # Equal Earth Asia-Pacific
  mutate(x=st_coordinates(.)[,1]/1e3,
         y=st_coordinates(.)[,2]/1e3,
         lrad=log(rad)) |>
  st_drop_geometry() |>
  arrange(track, date) |>
  group_by(track) |>
  mutate(x1=lag(x, d*1), y1=lag(y, d*1), lrad1=lag(lrad, d*1),
         x2=lag(x, d*2), y2=lag(y, d*2), lrad2=lag(lrad, d*2),
         x3=lag(x, d*3), y3=lag(y, d*3), lrad3=lag(lrad, d*3),
         x8=lag(x, d*8), y8=lag(y, d*8), lrad8=lag(lrad, d*8)) |>
  ungroup() |>
  rename(x0=x, y0=y) |>
  mutate(dx0=x0-x1, dy0=y0-y1,
         dx1=x1-x2, dy1=y1-y2,
         dx2=x2-x3, dy2=y2-y3,
         dx8=x1-x8, dy8=y1-y8,
         dist0=sqrt(dx0^2 + dy0^2),
         dist1=sqrt(dx1^2 + dy1^2),
         dist2=sqrt(dx2^2 + dy2^2),
         dist8=sqrt(dx8^2 + dy8^2),
         theta0=atan2(dy0, dx0),
         theta1=atan2(dy1, dx1),
         theta2=atan2(dy2, dx2),
         theta8=atan2(dy8, dx8)) %>%
  filter(complete.cases(.)) |>
  filter(dist0 < switch(as.character(d), "1"=150, "7"=300, "14"=500)) |>
  mutate(ldist0=log(dist0),
         ldist1=log(dist1),
         ldist2=log(dist2),
         ldist8=log(dist8)) |>
  mutate(altTheta0=if_else(theta0<0, theta0+pi, theta0-pi),
         altTheta1=if_else(theta1<0, theta1+pi, theta1-pi),
         altTheta2=if_else(theta2<0, theta2+pi, theta2-pi),
         altTheta8=if_else(theta8<0, theta8+pi, theta8-pi)) |>
  mutate(across(starts_with("ld"), ~c(scale(.x)))) 
ld_z <- tracks_df |> reframe(across(starts_with("dist"), ~c(log(mean(.x)), log(sd(.x)))))
saveRDS(ld_z, glue("data/cruise2024/aviso_ld_z_{d}d.rds"))
set.seed(1)
track_ids <- unique(tracks_df$track)
train_index <- sample(length(track_ids), length(track_ids)*0.8)
test_index <- seq_along(track_ids)[-train_index]
tracks_df <- tracks_df |>
  mutate(testSet=track %in% track_ids[test_index])
tracks_train <- tracks_df |> filter(!testSet)
tracks_test <- tracks_df |> filter(testSet)
 

ld_z <- readRDS(glue("data/cruise2024/aviso_ld_z_{d}d.rds"))
tracks_jc257 <- st_read("data/cruise2024/cruise2024_eddies.shp") %>%
  mutate(rad=as.numeric(sqrt(st_area(.)/pi))) |>
  st_centroid() |>
  mutate(track="jc257") |>
  add_lonlat() |>
  st_transform(8859) %>% # Equal Earth Asia-Pacific
  mutate(x=st_coordinates(.)[,1]/1e3,
         y=st_coordinates(.)[,2]/1e3,
         lrad=log(rad)) |>
  st_drop_geometry() |>
  arrange(track, date) |>
  group_by(track) |>
  mutate(x1=lag(x, d*1), y1=lag(y, d*1), lrad1=lag(lrad, d*1),
         x2=lag(x, d*2), y2=lag(y, d*2), lrad2=lag(lrad, d*2),
         x3=lag(x, d*3), y3=lag(y, d*3), lrad3=lag(lrad, d*3),
         x8=lag(x, d*8), y8=lag(y, d*8), lrad8=lag(lrad, d*8)) |>
  ungroup() |>
  rename(x0=x, y0=y) |>
  mutate(dx0=x0-x1, dy0=y0-y1,
         dx1=x1-x2, dy1=y1-y2,
         dx2=x2-x3, dy2=y2-y3,
         dx8=x1-x8, dy8=y1-y8,
         dist0=sqrt(dx0^2 + dy0^2),
         dist1=sqrt(dx1^2 + dy1^2),
         dist2=sqrt(dx2^2 + dy2^2),
         dist8=sqrt(dx8^2 + dy8^2),
         theta0=atan2(dy0, dx0),
         theta1=atan2(dy1, dx1),
         theta2=atan2(dy2, dx2),
         theta8=atan2(dy8, dx8)) %>%
  filter(complete.cases(.)) |>
  mutate(ldist0=log(dist0),
         ldist1=log(dist1),
         ldist2=log(dist2),
         ldist8=log(dist8)) |>
  mutate(altTheta0=if_else(theta0<0, theta0+pi, theta0-pi),
         altTheta1=if_else(theta1<0, theta1+pi, theta1-pi),
         altTheta2=if_else(theta2<0, theta2+pi, theta2-pi),
         altTheta8=if_else(theta8<0, theta8+pi, theta8-pi)) |>
  mutate(ldist0=(ldist0-ld_z$dist0[1])/ld_z$dist0[2],
         ldist1=(ldist1-ld_z$dist1[1])/ld_z$dist1[2],
         ldist2=(ldist2-ld_z$dist2[1])/ld_z$dist2[2],
         ldist8=(ldist8-ld_z$dist8[1])/ld_z$dist8[2]) 
  tracks_train <- bind_rows(tracks_train, tracks_jc257)


# out <- brm(bf(ldist0 ~ s(lon1,lat1) + ldist1 + (1|track)),
#            data=tracks_train, cores=4, chains=4, family=skew_normal(),
#            prior=prior(normal(0, 1), class="b"))
# summary(out)
# pp_check(out)
# conditional_effects(out)
# conditional_smooths(out, )
# bayes_R2(out)
# tracks_train |>
#   mutate(pred=colMeans(posterior_epred(out))) |>
#   ggplot(aes(pred, ldist0)) + geom_point()
# 
# out <- brm(bf(altTheta0 ~ s(lon1,lat1) + altTheta1 + (1|track)),
#            data=tracks_train, chains=3, cores=3, family=von_mises(),
#            prior=c(prior(normal(0, 3), class="Intercept", lb=-3.141593, ub=3.141593),
#                    prior(normal(0, 1), class="b"),
#                    prior(normal(0, 1), class="kappa", lb=0)))
# summary(out)
# pp_check(out)
# conditional_effects(out)
# conditional_smooths(out, surface=T)
# bayes_R2(out)
# tracks_train |>
#   mutate(pred=colMeans(posterior_epred(out))) |>
#   ggplot(aes(pred, altTheta0)) + geom_point()
# 
# 
# out <- brm(bf(lrad ~ lon1 + lat1 + lrad1 + lrad2 + (1|track)), 
#            data=tracks_train, cores=4, chains=4,
#            prior=prior(normal(0, 1), class="b"))
# summary(out)
# pp_check(out)
# conditional_effects(out)  
# bayes_R2(out)
# tracks_train |>
#   mutate(pred=colMeans(posterior_epred(out))) |>
#   ggplot(aes(pred, lrad)) + geom_point()


bf_ldist0 <- bf(ldist0 ~ s(x1,y1) + ldist1 + ldist2 + ldist8 + 
                  (1 + ldist1 + ldist2 + ldist8 | i | track)) + skew_normal()
bf_theta0 <- bf(altTheta0 ~ s(x1,y1) + altTheta1 + altTheta2 + altTheta8 + 
                  (1 + altTheta1 + altTheta2 + altTheta8 | j | track)) + von_mises()
bf_lrad <- bf(lrad ~ s(x1,y1) + lrad1 + lrad2 + lrad8 + 
                (1 + lrad1 + lrad2 + lrad8 | k | track)) + gaussian()
out_mv <- brm(bf_ldist0 + bf_theta0 + bf_lrad, 
              data=tracks_train, cores=3, chains=3, threads=threading(20),
              refresh=10, 
              control=list(max_treedepth=15),
              prior=c(prior(normal(0, 1), class="b", resp="altTheta0"),
                      prior(normal(0, 1), class="b", resp="ldist0"),
                      prior(normal(0, 1), class="b", resp="lrad"),
                      prior(normal(0, 1), class="sd", resp="ldist0", group="track"),
                      prior(normal(0, 1), class="sd", resp="altTheta0", group="track"),
                      prior(normal(0, 1), class="sd", resp="lrad", group="track"),
                      prior(normal(0, 1), class="Intercept", resp="ldist0"),
                      prior(normal(12, 2), class="Intercept", resp="lrad"),
                      prior(normal(0, 3), class="Intercept", resp="altTheta0", lb=-3.141593, ub=3.141593),
                      prior(normal(0, 1), class="kappa", resp="altTheta0", lb=0)))
saveRDS(out_mv, glue("data/cruise2024/aviso_{d}day-lag_brm.rds"))








out_mv <- readRDS("data/cruise2024/aviso_7day-lag_brm.rds")

out_mv
pp_check(out_mv, resp="ldist0")
pp_check(out_mv, resp="altTheta0")
pp_check(out_mv, resp="lrad")

conditional_effects(out_mv)
plot(conditional_smooths(out_mv, plot=F), stype="raster")

bayes_R2(out_mv)


fits <- posterior_epred(out_mv)
fit_df <- tracks_train |>
  mutate(p_ldist0=colMeans(fits[,,1]),
         p_altTheta0=colMeans(fits[,,2]),
         p_lrad=colMeans(fits[,,3])) |>
  mutate(p_theta0=p_altTheta0 + pi*if_else(p_altTheta0>0, -1, 1),
         p_dist0=pmin(exp(p_ldist0*ld_z$dist0[2] + ld_z$dist0[1]), 5e2),
         p_dx0=p_dist0 * cos(p_theta0) * 1e3,
         p_dy0=p_dist0 * sin(p_theta0) * 1e3,
         p_x0=x1+p_dx0,
         p_y0=y1+p_dy0)

ggplot(fit_df) +
  geom_path(aes(x0, y0, group=track), linewidth=1) +
  geom_path(aes(p_x0, p_y0, group=track), colour="red")

fit_df |>
  ggplot(aes(p_ldist0, ldist0)) + geom_point()
fit_df |>
  ggplot(aes(p_dist0, dist0)) + geom_point()
fit_df |>
  ggplot(aes(p_dx0, dx0)) + geom_point()




preds <- posterior_epred(out_mv, newdata=tracks_test, allow_new_levels=T)
pred_df <- tracks_test |>
  mutate(p_ldist0=colMeans(preds[,,1]),
         p_altTheta0=colMeans(preds[,,2]),
         p_lrad=colMeans(preds[,,3])) |>
  mutate(p_theta0=p_altTheta0 + pi*if_else(p_altTheta0>0, -1, 1),
         p_dist0=pmin(exp(p_ldist0*ld_z$dist0[2] + ld_z$dist0[1]), 5e2),
         p_dx0=p_dist0 * cos(p_theta0) * 1e3,
         p_dy0=p_dist0 * sin(p_theta0) * 1e3,
         p_x0=x1+p_dx0,
         p_y0=y1+p_dy0,
         p_rad=exp(p_lrad))

ggplot(pred_df) +
  geom_point(aes(x0, y0, group=track), size=2) +
  geom_point(aes(p_x0, p_y0, group=track), colour="red") +
  facet_wrap(~track, scales="free")

pred_df |>
  ggplot(aes(p_x0, x0)) + geom_point() + geom_abline()
pred_df |>
  ggplot(aes(p_y0, y0)) + geom_point() + geom_abline()

pred_df |>
  ggplot(aes(p_ldist0, ldist0)) + geom_point() + geom_abline()
pred_df |>
  ggplot(aes(p_dist0, dist0)) + geom_point() + geom_abline()

pred_df |>
  ggplot(aes(p_lrad, lrad)) + geom_point() + geom_abline()
pred_df |>
  ggplot(aes(p_rad, rad)) + geom_point() + geom_abline()
pred_df |>
  ggplot(aes(p_altTheta0, altTheta0)) + geom_point() + geom_abline()

pred_df |>
  ggplot(aes((p_ldist0 - ldist0)/abs(ldist0))) + geom_histogram()


