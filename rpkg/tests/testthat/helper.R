london <- function() {
  sf::st_read(
    system.file("extdata", "london.geojson", package = "gridmapper"),
    quiet = TRUE
  )
}

js_reference <- function() {
  jsonlite::read_json(test_path("fixtures", "js-reference.json"))
}

as_num_matrix <- function(x) {
  do.call(rbind, lapply(x, \(r) as.numeric(unlist(r))))
}
