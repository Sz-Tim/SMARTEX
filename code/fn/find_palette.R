



find_palette <- function(var, lims, type="fill") {
  
  if(type != "fill") {
    stop("Only fill is coded so far -- add desired palette")
  }
  
  if(var=="sla") {
    pal <- scale_fill_continuous_divergingx(
      name='SSH\nanomaly (m)',
      palette='RdBu', mid=0, limits=lims, l3 = 0, p3 = .8, p4 = .6, rev=T
    ) 
  }
  
  if(var %in% c("EKE", "EKE_std")) {
    pal <- scale_fill_viridis_c(
      ifelse(var=="EKE",
             expression(paste("EKE (cm" ^2, "s" ^-2, ")")),
             expression(z['dep'](EKE))), 
      option="inferno", limits=lims
      )
  }
  
  if(var %in% c("u", "u_std")) {
    pal <- scale_fill_continuous_divergingx(
      name=ifelse(var=="u", "u (m/s)", expression(z['dep'](u))), 
      palette="PuOr", mid=0, rev=T, limits=lims
      )
  }
  
  if(var %in% c("v", "v_std")) {
    pal <- scale_fill_continuous_divergingx(
      name=ifelse(var=="v", "v (m/s)", expression(z['dep'](v))), 
      palette="RdBu", mid=0, rev=T, limits=lims)
  }
  
  if(var %in% c("u2v2", "u2v2_std")) {
    pal <- scale_fill_viridis_c(
      name=ifelse(var=="u2v2", 
                  expression(paste(u^2 + v^2, "(m" ^2, "s" ^-2, ")")), 
                  expression(z['dep'](u^2+v^2))), 
      option="inferno", limits=lims)
  }
  
  if(var %in% c("uv", "uv_std")) {
    pal <- scale_fill_viridis_c(
      name=ifelse(var=="uv", 
                  "Current\nspeed (m/s)", 
                  expression(z['dep'](sqrt(u^2+v^2)))), 
      option="inferno", limits=lims, oob=scales::squish)
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