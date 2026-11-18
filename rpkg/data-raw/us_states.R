# US state centroids, from demo/usStates.json of the JavaScript library.
states <- jsonlite::fromJSON("../demo/usStates.json")
idx <- match(states$name, datasets::state.name)
us_states <- data.frame(
  name = states$name,
  abbr = datasets::state.abb[idx],
  region = as.character(datasets::state.region[idx]),
  lon = states$lon,
  lat = states$lat
)
usethis::use_data(us_states, overwrite = TRUE)
