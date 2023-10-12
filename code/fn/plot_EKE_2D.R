# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


plot_EKE_2D <- function(cmems_df, bathy, ccz, POI, fill_var, fill_title,
                         xlim=c(NA_real_, NA_real_), ylim=c(NA_real_, NA_real_), 
                         out_f=NA, out_dim=c(11, 5)) {
  
  p <- cmems_df |> 
    group_by(lat, lon) |> 
    summarise(EKE_mn=mean(EKE), 
              EKE_var=var(EKE),
              EKE_sd=sd(EKE)) |> 
    ungroup() |>
    ggplot() + 
    geom_sf(data=bathy, aes(fill=elev), colour=NA, alpha=0.5) +
    scale_fill_distiller(type="seq", palette="Greys", direction=1, guide="none") +
    new_scale_fill() +
    geom_tile(aes(lon, lat, fill=.data[[fill_var]])) + 
    scale_fill_viridis_c(fill_title, 
                         option="inferno") + 
    geom_sf(data=ccz, linewidth=0.25, colour="grey") +
    geom_sf(data=ccz |> filter(Contractor=="UKSRL"), linewidth=0.75, colour="grey") +
    geom_point(data=POI, aes(lon, lat, shape=name), colour="darkslategray1", size=3) +
    scale_shape_manual("", values=c(5,19,3,4)) +
    xlim(xlim[[1]], xlim[[2]]) + ylim(ylim[[1]], ylim[[2]]) +
    ggtitle(paste(min(cmems_df$date), "to", max(cmems_df$date))) +
    guides(fill=guide_colourbar(order=1), 
           shape=guide_legend(order=2, override.aes=list(colour="black"))) +
    theme(axis.title=element_blank(), 
          legend.position="bottom", 
          legend.key.height=unit(0.25, "cm"), 
          legend.key.width=unit(1.25, "cm"))
  if(is.na(out_f)) {
    return(p)
  } else {
    ggsave(out_f, p, width=out_dim[1], height=out_dim[2], units="in", dpi=200)
  }
}