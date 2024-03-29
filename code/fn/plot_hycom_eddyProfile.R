


plot_hycom_eddyProfile <- function(df, depthLim=4000, lon_diff, POI, 
                                   fill_var, fill_lim, alpha_lim=c(0,1), title, 
                                   xlim=c(NA_real_, NA_real_), ylim=c(NA_real_, NA_real_), 
                                   out_f=NA, out_dim=NULL, tracks=NULL,
                                   darkLines=TRUE) {
  
  fill_pal <- find_palette(fill_var, fill_lim, "fill")
  if(darkLines) {
    cols <- list(ccz="grey30",
                 POI="black",
                 active="black",
                 inactive="grey30")
  } else {
    cols <- list(ccz="grey",
                 POI="darkslategray3",
                 active="grey95",
                 inactive="grey80")
  }
  
  p <- df |>
    filter(depth < depthLim) |>
    ggplot() +
    geom_line(data=tracks, aes(lon, y=0), linewidth=2, colour="grey40") +
    {
      if(fill_var != "uvDir") {
        geom_rect(aes(xmin=lon-lon_diff, xmax=lon+lon_diff, 
                      ymin=depth0, ymax=depth1,
                      fill=.data[[fill_var]]), colour=NA)
      } else {
        geom_rect(aes(xmin=lon-lon_diff, xmax=lon+lon_diff, 
                      ymin=depth0, ymax=depth1,
                      fill=.data[[fill_var]], alpha=EKE_std), colour=NA)
      }
    } + 
    geom_point(data=POI, aes(lon, y=0), shape=4, size=3, colour=cols$POI) +
    scale_alpha_continuous(expression(z['dep'](EKE)), range=c(0, 1), limits=alpha_lim) +
    fill_pal + 
    scale_y_reverse() +
    labs(x="Longitude", y="Depth (m)", title=title) +
    guides(fill=guide_colourbar(order=1), 
           shape=guide_legend(order=2, override.aes=list(colour="black"), nrow=1)) +
    theme(legend.position="bottom")
  if(is.na(out_f)) {
    return(p)
  } else {
    ggsave(out_f, p, width=out_dim[1], height=out_dim[2], units="in", dpi=200)
  }
}
