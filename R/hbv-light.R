
#' Title
#'
#' @param hbv_light_dir
#'
#' @return
#' @export
#'
#' @examples
parse_hbv_light <- function(hbv_light_dir) {
  # not supported: sub-daily values, 360day calendar
  # evap 12-months
  ptq_file <- list.files(hbv_light_dir,full.names=TRUE,pattern="(?i)ptq\\.txt$")[1]
  evap_file <- list.files(hbv_light_dir,full.names=TRUE,pattern="(?i)evap\\.txt$")[1]
  clarea_file <- list.files(hbv_light_dir,full.names=TRUE,pattern="(?i)clarea\\.xml$")[1]
  ptq <- read_ptq(ptq_file)
  ep <- read_evap(evap_file,ref_ts = ptq)
  clarea <- read_clarea(clarea_file)
  return(list(prec=ptq[,1], airt=ptq[,2], ep=ep, obs=ptq[,3], area=clarea$area_elev, elev_zones=clarea$zones))
}


#' Extract settings from HBV-Light's CLAREA file
#'
#' Reads elevation zones and their relative area into vectors. Caution: Usage of subcatchments is currently not supported. Contatc the author if you require this functionality.
#' @param clarea_file File path to CLAREA.xml
#'
#' @return List containing the following named items:
#' \item{zones}{Numeric vector of the elevation zones (m.a.s.l)}
#' \item{area_elev}{Numeric vector of the eleavtion zone's relative}
#' \item{area_vegetation}{Matrix describing the relative area of the vegetation zones (vegetation zones in columns))}
#' @export
#'
#' @examples
#' @author Joschka Thurner  \email{joschka.thurner@th-koeln.de}
#' @importFrom magrittr %>%
#' @import xml2
read_clarea <- function(clarea_file) {
  clarea <- xml2::read_xml(clarea_file)
  zones <- clarea %>%
    xml2::xml_find_first("ElevationZoneHeight") %>%
    xml2::xml_children() %>%
    xml2::xml_double()
  area_vegetation <- clarea %>%
    xml2::xml_find_all(".//VegetationZone") %>%
    lapply(FUN=xml2::xml_children) %>%
    sapply(FUN=xml2::xml_double)

  return(list("zones"=zones,"area_elev"=rowSums(area_vegetation),"area_vegetation"=area_vegetation))
}

#' Title
#'
#' @param ptq_file
#'
#' @return
#' @export
#'
#' @examples
read_ptq <- function(ptq_file) {

  ptq <- as.zoo(read.delim.zoo(ptq_file, skip = 2, header = FALSE, col.names = c("date","prec","airt","qobs"),
                               format = "%Y%m%d",tz = NULL))
  # hbv-light codes missing qobs as negative values
  ptq$qobs[ptq$qobs < 0] <- NA
  return(ptq)
  # ptq <- readr::read_delim(ptq_file,
  #            "\t", escape_double = FALSE, col_types = readr::cols(date = readr::col_date(format = "%Y%m%d")),
  #            comment = "#", trim_ws = TRUE, skip = 1)
  # # return(ptq)
  # # return(dplyr::mutate_all(ptq,dplyr::funs(replace(., .<0, NA))))
  # return(mutate(ptq,qobs=replace(qobs,which(qobs<0),NA)))
}


#' Title
#'
#'no subcatchemnt support
#'
#' @param evap_file
#' @param ref_ts
#'
#' @return
#' @export
#'
#' @examples
read_evap <- function(evap_file, ref_ts=NULL) {
  # TODO: handle yearly daily mean and daily values (rep/add index)
  ep <- read.delim(evap_file)[, 1]
  if (!is.null(ref_ts)) {
    ep <- monthly2daily(start(ref_ts), end(ref_ts), ep, FUN = na.spline)
  }
  return(ep)
}




#' #' Extract settings from HBV-Light's SIMULATION file
#' #'
#' #' Reads elevation zones and their relative area into vectors. Caution: Usage of subcatchments is currently not supported. Contatc the author if you require this functionality.
#' #' @param clarea_file File path to Simulation.xml
#' #'
#' #' @return List containing the vectors "zones" and "area", corresponding to the elevation zones and their relative area as set in HBV-Light.
#' #' @export
#' #'
#' #' @examples
#' #' @author Joschka Thurner  \email{joschka.thurner@th-koeln.de}
#' read_simulation <- function(simulation_file) {
#'
#' }
#'
#'
#' #' Title
#' #'
#' #' @param parameter_file
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' read_parameter <- function(parameter_file) {
#'   hbv_pars <- xml2::read_xml(parameter_file)
#'   par_type <- xml2::xml_name(hbv_pars)
#'   if (par_type == "Catchment") {
#'     catchment_pars <- xml2::xml_find_first(hbv_pars,".//SubCatchmentParameters") %>%
#'       xml2::xml_children() %>%
#'       xml2::xml_double(.)
#'
#'       setNames(object=xml2::xml_double(.),nm=xml2::xml_name(.))
#'     # catchment_pars <- setNames(xml_double(xml_children(d)),xml_name(xml_children(d)))
#'
#'   }
#'
#' }
#'


