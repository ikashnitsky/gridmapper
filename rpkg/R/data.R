#' Locations of the US states
#'
#' Representative point locations of the 50 US states.
#'
#' @format A data frame with 50 rows and 5 columns:
#' \describe{
#'   \item{name}{State name.}
#'   \item{abbr}{Two-letter postal abbreviation.}
#'   \item{region}{Census region, from [datasets::state.region].}
#'   \item{lon, lat}{Longitude and latitude in degrees.}
#' }
#' @source The demo data of the 'gridmapper' JavaScript library.
"us_states"
