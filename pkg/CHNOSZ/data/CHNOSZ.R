# CHNOSZ/data/CHNOSZ.R
# clear system settings (basis/species)

# we only work if the CHNOSZ environment exists
if(!"CHNOSZ" %in% search()) {
  message("data(CHNOSZ): please run data(thermo) first")
} else {

  local({
    # get thermo from CHNOSZ environment
    thermo <- get("thermo", "CHNOSZ")
    # set basis,species components
    thermo$basis <- NULL
    thermo$species <- NULL
    # place it in CHNOSZ environment
    assign("thermo", thermo, "CHNOSZ")
  })

  message("data(CHNOSZ): cleared basis and species settings")

}
