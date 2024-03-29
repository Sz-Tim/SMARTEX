

plot_hycom_depthPanels <- function(df, depths, bathy=NULL, ccz, POI, 
                                   fill_var, fill_lim, alpha_lim=c(0,1), title, 
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
    filter(depth %in% depths) |>
    ggplot() + 
    {if(!is.null(bathy)) {
      geom_sf(data=bathy, aes(fill=elev), colour=NA, alpha=0.5)
    }} +
    {if(!is.null(bathy)) {
      scale_fill_distiller(type="seq", palette="Greys", direction=1, guide="none")
    }} +
    {if(!is.null(bathy)) {
      new_scale_fill()
    }} +
    {if(fill_var != "uvDir") {
      geom_raster(aes(lon, lat, fill=.data[[fill_var]]), 
                alpha=ifelse(!is.null(bathy), 0.5, 1))
    } else {
      geom_raster(aes(lon, lat, fill=.data[[fill_var]], alpha=EKE_std))
    }} + 
    scale_alpha_continuous(expression(z['dep'](EKE)), range=c(0.5, 1), limits=alpha_lim) +
    fill_pal + 
    {if(!is.null(tracks)) {
      geom_path(data=tracks, aes(lon, lat, group=track, colour=active, linewidth=active))
    }} +
    scale_colour_manual(values=c("FALSE"=cols$inactive, "TRUE"=cols$active), guide="none") +
    geom_sf(data=ccz, linewidth=0.25, colour=cols$ccz) +
    geom_point(data=POI, aes(lon, lat, shape=name), colour=cols$POI) +
    scale_shape_manual("", values=c(5,19,3,4)) +
    scale_linewidth_manual("", values=c("FALSE"=0.4, "TRUE"=0.8), guide="none") + 
    labs(title=title) + 
    xlim(xlim[[1]], xlim[[2]]) + ylim(ylim[[1]], ylim[[2]]) +
    facet_grid(depth~.) +
    guides(fill=guide_colourbar(order=1), 
           shape=guide_legend(order=2, override.aes=list(colour="black"), ncol=1)) +
    theme(axis.title=element_blank())
  if(is.na(out_f)) {
    return(p)
  } else {
    ggsave(out_f, p, width=out_dim[1], height=out_dim[2], units="in", dpi=200)
  }
  
}
