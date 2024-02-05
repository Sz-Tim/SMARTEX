



find_palette <- function(var, lims, type="fill") {
  
  library(colorspace)
  
  if(type != "fill") {
    if(var == "uvDir") {
      pal <- scale_colour_gradientn(
        name="Current\ndirection",
        colours=readRDS("data/cmr_cmaps.RDS")$infinity, limits=c(-pi, pi),
        breaks=c(-pi, -pi/2, 0, pi/2, pi),
        labels=c("W", "S", "E", "N", "W"))
    }
  }
  
  if(var %in% c("sla", "el_diff")) {
    pal <- scale_fill_continuous_divergingx(
      name=if_else(var=="sla", 'SSH\nanomaly (m)', 'Difference\nHY-CM (m)'),
      palette='RdBu', mid=0, limits=lims, l3 = 0, p3 = .8, p4 = .6, rev=T
    ) 
  }
  
  if(var=="surf_el") {
    pal <- scale_fill_viridis_c(
      name='SSH (m)', limits=lims
    ) 
  }
  
  if(var %in% c("EKE", "EKE_std", "EKE_mn")) {
    pal <- scale_fill_viridis_c(
      ifelse(var %in% c("EKE", "EKE_mn"),
             expression(paste("EKE (cm" ^2, "s" ^-2, ")")),
             expression(z['dep'](EKE))), 
      option="turbo", begin=0.05, end=0.95, limits=lims
      )
  }
  
  if(var %in% c("EKE_sd")) {
    pal <- scale_fill_viridis_c(
      "Interannual sd\nin EKE", 
      option="cividis", limits=lims
    )
  }
  
  if(var %in% c("EKE_cv", "uv_cv")) {
    pal <- scale_fill_viridis_c(
      ifelse(var == "EKE_cv",
             "Interannual CV\nin EKE",
             "Interannual CV\nin current speed"), 
      option="cividis", limits=lims
    )
  }
  
  if(var %in% c("u", "u_std", "ugos", "u_diff")) {
    pal <- scale_fill_continuous_divergingx(
      name=switch(var,
                  "u"="u (m/s)",
                  "ugos"="ugos (m/s)",
                  "u_std"=expression(z['dep'](u)),
                  "u_diff"="Difference\nHY-CM (m/s)"),
      palette="PuOr", mid=0, rev=T, limits=lims
      )
  }
  
  if(var %in% c("v", "v_std", "vgos", "v_diff")) {
    pal <- scale_fill_continuous_divergingx(
      name=switch(var,
                  "v"="v (m/s)",
                  "vgos"="vgos (m/s)",
                  "v_std"=expression(z['dep'](v)),
                  "v_diff"="Difference\nHY-CM (m/s)"),
      palette="RdBu", mid=0, rev=T, limits=lims)
  }
  
  if(var %in% c("u2v2", "u2v2_std")) {
    pal <- scale_fill_viridis_c(
      name=ifelse(var=="u2v2", 
                  expression(paste(u^2 + v^2, "(m" ^2, "s" ^-2, ")")), 
                  expression(z['dep'](u^2+v^2))), 
      option="turbo", begin=0.05, end=0.95, limits=lims)
  }
  
  if(var %in% c("uv", "uv_std", "uv_mn")) {
    pal <- scale_fill_viridis_c(
      name=switch(var,
                  "uv"="Current\nspeed (m/s)",
                  "uv_mn"="Current\nspeed (m/s)",
                  "uv_std"=expression(z['dep'](sqrt(u^2+v^2)))),
      option="turbo", begin=0.05, end=0.95, limits=lims, oob=scales::squish)
  }
  
  if(var == "uv_diff") {
    pal <- scale_fill_continuous_divergingx(
      name="Difference\nHY-CM (m/s)",
      palette="PRGn", mid=0, rev=T, limits=lims)
  }
  
  if(var %in% c("uv_sd")) {
    pal <- scale_fill_viridis_c(
      "Interannual sd\nin current speed", 
      option="cividis", limits=lims
    )
  }
  
  if(var == "uvDir") {
    pal <- scale_fill_gradientn(
      name="Current\ndirection",
      colours=readRDS("data/cmr_cmaps.RDS")$infinity, limits=c(-pi, pi),
      breaks=c(-pi, -pi/2, 0, pi/2, pi),
      labels=c("W", "S", "E", "N", "W"))
  }
  
  return(pal)
}