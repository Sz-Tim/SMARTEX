# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


plot_cmems <- function(df, bathy, ccz, POI, fill_var, fill_lim, title, 
                       xlim=c(NA_real_, NA_real_), ylim=c(NA_real_, NA_real_), 
                       out_f=NA, out_dim=c(),
                       tracks=NULL) {
  
  if(fill_var=="sla") {
    fill_pal <- scale_fill_continuous_divergingx(
      name='SSH\nanomaly (m)',
      palette='RdBu', mid=0, limits=fill_lim, l3 = 0, p3 = .8, p4 = .6, rev=T
    )
  }
  if(fill_var=="EKE") {
    fill_pal <- scale_fill_viridis_c(expression(paste("EKE (cm" ^2, "s" ^-2, ")")), 
                                     option="inferno", limits=fill_lim)
  }
  
  p <- ggplot(df) + 
    geom_sf(data=bathy, aes(fill=elev), colour=NA, alpha=0.5) +
    scale_fill_distiller(type="seq", palette="Greys", direction=1, guide="none") +
    new_scale_fill() +
    geom_tile(aes(lon, lat, fill=.data[[fill_var]]), alpha=0.5, colour=NA) + 
    fill_pal + 
    {if(!is.null(tracks)) {
      geom_path(data=tracks, aes(lon, lat, group=track, colour=active, linewidth=active))
    }} +
    scale_colour_manual(values=c("FALSE"="grey", "TRUE"="black"), guide="none") +
    geom_sf(data=ccz, linewidth=0.25, colour="grey30") +
    geom_point(data=POI, aes(lon, lat, shape=name)) +
    scale_shape_manual("", values=c(1,4)) +
    scale_linewidth_manual("", values=c("FALSE"=0.4, "TRUE"=0.8), guide="none") + 
    labs(title=title) + 
    xlim(xlim[[1]], xlim[[2]]) + ylim(ylim[[1]], ylim[[2]]) +
    theme(axis.title=element_blank())
  if(is.na(out_f)) {
    return(p)
  } else {
    ggsave(out_f, p, width=out_dim[1], height=out_dim[2], units="in", dpi=200)
  }
}
