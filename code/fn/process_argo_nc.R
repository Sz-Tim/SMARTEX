#' Process ARGOS nc files
#'
#' @param f_nc filename
#'
#' @return
#' @export
#'
#' @examples
process_argo_nc <- function(f_nc) {
  f_parts <- str_split(f_nc, "/")[[1]]
  i_nc <- nc_open(f_nc) 
  argos_df <- tibble(
    PN=c(ncvar_get(i_nc, "PLATFORM_NUMBER")),
    CN=c(ncvar_get(i_nc, "CYCLE_NUMBER")),
    FSN=c(ncvar_get(i_nc, "FLOAT_SERIAL_NO")),
    juld=c(ncvar_get(i_nc, "JULD")),
    mtime=ymd("1950-01-01") + juld,
    lon=c(ncvar_get(i_nc, "LONGITUDE")),
    lat=c(ncvar_get(i_nc, "LATITUDE")),
    fnm=f_nc,
    dir_id=str_split_fixed(last(f_parts), "_", 2)[,1]
  )
  if(nrow(argos_df)>1) {
    argos_df <- argos_df %>%
      filter(!duplicated(.))
  }
  n_obs <- i_nc$dim$N_LEVELS$len
  n_prof <- i_nc$dim$N_PROF$len
  prof_lu <- tibble(df_name=c("Pi", "Ti", "Si", 
                              "P", "T", "S"),
                    nc_name=c("PRES", "TEMP", "PSAL", 
                              "PRES_ADJUSTED", "TEMP_ADJUSTED", "PSAL_ADJUSTED"))
  prof_ls <- map(prof_lu$df_name, NA_real_) |> set_names(prof_lu$df_name)
  for(i in seq_along(prof_ls)) {
    dat_i <- NULL
    try({
      dat_i <- ncvar_get(i_nc, prof_lu$nc_name[i]) 
    })
    if(!is.null(dat_i)) {
      prof_ls[[i]] <- tibble(val=c(dat_i), 
                             obs=rep(1:n_obs, n_prof),
                             profile=rep(1:n_prof, each=n_obs)) |>
        rename_with(~prof_lu$df_name[i], val)
    }
  }
  nc_close(i_nc)
  
  prof_df <- bind_cols(argos_df, reduce(prof_ls, full_join, by=c("obs", "profile")))
  return(prof_df)
}