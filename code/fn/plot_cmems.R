# 
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


plot_cmems <- function(df, bathy=NULL, land=NULL, ccz, POI, 
                       fill_var, fill_lim, alpha_lim, title, 
                       xlim=c(NA_real_, NA_real_), ylim=c(NA_real_, NA_real_), 
                       out_f=NA, out_dim=NULL,
                       tracks=NULL, darkLines=TRUE) {
  
  fill_pal <- find_palette(fill_var, fill_lim, "fill")
  if(darkLines) {
    cols <- list(ccz="grey30",
                 POI="black",
                 active="black",
                 inactive="grey30")
  } else {
    cols <- list(ccz="grey",
                 POI="darkslategray1",
                 active="grey95",
                 inactive="grey80")
  }
  
  p <- df |>
    ggplot() + 
    {if(!is.null(bathy)) {
      geom_stars(data=bathy, downsample=10) }
      # geom_sf(data=bathy, aes(fill=elev), colour=NA, alpha=0.5) }
    } +
    {if(!is.null(bathy)) {
      scale_fill_continuous_divergingx(h1=0, h2=NA, h3=0, 
                                       c1=0, c2=NA, c3=0, 
                                       l1=0, l2=90, l3=0, 
                                       p1=1.3, p2=NA, p3=1.3, p4=0.75, 
                                       cmax1=NA, cmax2=NA, guide="none") }
      # scale_fill_distiller(type="seq", palette="Greys", direction=1, guide="none")   } 
    } +
    {if(!is.null(bathy)) {
      new_scale_fill()}
    } +
    {if(!is.null(land)) {
      geom_sf(data=land) } 
    } +
    # geom_tile(aes(lon, lat, fill=.data[[fill_var]]), alpha=1, colour=NA) +
    geom_tile(aes(lon, lat, fill=.data[[fill_var]], alpha=sqrt(abs(sla))), colour=NA) +
    scale_alpha_continuous(limits=alpha_lim, range=c(0.5, 1), guide="none") +
    fill_pal + 
    {if(!is.null(tracks)) {
      geom_path(data=tracks, aes(lon, lat, group=track, colour=active, linewidth=active))
    }} +
    scale_colour_manual(values=c("FALSE"=cols$inactive, "TRUE"=cols$active), guide="none") +
    geom_sf(data=ccz, linewidth=0.25, colour=cols$ccz) +
    geom_sf(data=ccz |> filter(Contractor=="UKSRL"), linewidth=0.75, colour=cols$ccz) +
    geom_point(data=POI, aes(lon, lat, shape=name), colour=cols$POI, size=3) +
    scale_linewidth_manual("", values=c("FALSE"=0.4, "TRUE"=0.8), guide="none") + 
    scale_shape_manual("", values=c(5,19,3,4)) +
    labs(title=title) + 
    xlim(xlim[[1]], xlim[[2]]) + ylim(ylim[[1]], ylim[[2]]) +
    guides(fill=guide_colourbar(order=1), 
           shape=guide_legend(order=2, override.aes=list(colour="black"), nrow=2)) +
    theme(axis.title=element_blank(), 
          legend.position="bottom", 
          legend.key.height=unit(0.4, "cm"), 
          legend.key.width=unit(1.25, "cm"))
  if(is.na(out_f)) {
    return(p)
  } else {
    ggsave(out_f, p, width=out_dim[1], height=out_dim[2], units="in", dpi=200)
  }
}
