## ---- include = FALSE---------------------------------------------------------
env_present <- slendr:::is_slendr_env_present()

knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  dpi = 60,
  eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS"
)

## -----------------------------------------------------------------------------
library(slendr)

init_env()

## -----------------------------------------------------------------------------
map <- world(
  xrange = c(0, 1000),
  yrange = c(0, 1000),
  landscape = "blank"
)

## -----------------------------------------------------------------------------
create_pop <- function(i, n_side, map, N, radius) {
  # get dimensions of the world map
  dim <- c(diff(attr(map, "xrange")), diff(attr(map, "yrange")))

  # position of the i-th population on the two-dimensional lattice grid
  coords <- c((i - 1) %% n_side, (i - 1) %/% n_side)
  center <- coords / n_side * dim + dim / (2 * n_side)

  pop <- tryCatch({
    population(
      name = sprintf("pop%d", i),
      N = N,
      time = 1,
      map = map,
      center = center + c(attr(map, "xrange")[1], attr(map, "yrange")[1]),
      radius = radius
    )
  }, error = function(e) NULL)

  pop
}

## -----------------------------------------------------------------------------
n <- 5

populations <-
  seq(1, n * n) %>%
  lapply(create_pop, n_side = n, map = map, N = 100, radius = 40)

## ----demes_grid---------------------------------------------------------------
do.call(plot_map, populations) + ggplot2::theme(legend.position = "none")

## -----------------------------------------------------------------------------
set_geneflow <- function(i, n_side, rate, start, end, populations) {
  pop <- populations[[i]]

  # get the position of the i-th population on the n*n grid
  coords <- c((i - 1) %% n_side, (i - 1) %/% n_side)

   # get coordinates of the i-th population's neighbors on the grid
  neighbor_pos <- list(
    c(coords[1] - 1, coords[2]), c(coords[1] + 1, coords[2]),
    c(coords[1], coords[2] + 1), c(coords[1], coords[2] - 1)
  )

  # generate geneflow events for population coordinates inside the grid
  geneflows <- lapply(neighbor_pos, function(pos) {
    if (any(pos < 0 | pos >= n_side)) return(NULL)
    neighbor <- populations[[pos[2] * n_side + pos[1] + 1]]
    if (is.null(neighbor)) return(NULL)

    rbind(
      gene_flow(from = pop, to = neighbor, rate = rate, start = start, end = end, overlap = FALSE),
      gene_flow(from = neighbor, to = pop, rate = rate, start = start, end = end, overlap = FALSE)
    )
  }) %>%
    do.call(rbind, .)

  geneflows
}

## -----------------------------------------------------------------------------
set_geneflow(1, n, rate = 0.1, start = 2, end = 1000, populations)

## -----------------------------------------------------------------------------
geneflows <-
  seq(1, n * n) %>%
  lapply(set_geneflow, n, rate = 0.05, start = 2, end = 1000, populations) %>%
  do.call(rbind, .) %>%
  unique # filter out duplicate events due to symmetries

## -----------------------------------------------------------------------------
nrow(geneflows)

## -----------------------------------------------------------------------------
model <- compile_model(
  populations = populations, gene_flow = geneflows,
  generation_time = 1, resolution = 10,
  competition = 10, mating = 10, dispersal = 10,
  simulation_length = 1000
)

## -----------------------------------------------------------------------------
ts <- slim(model, sequence_length = 10000, recombination_rate = 0) # simulate a single 10kb locus
ts

## ---- results = FALSE, warning = FALSE----------------------------------------
map <- world(
  xrange = c(-25, 55),
  yrange = c(-32, 35),
  crs = 4326
)

n <- 20

populations <-
  seq(1, n * n) %>%
  lapply(create_pop, n_side = n, map = map, N = 100, radius = 1.5)

## -----------------------------------------------------------------------------
continent <- region(
  map = map, polygon = list(
    c(-10, 35), c(-20, 20), c(-15, 8), c(-10, 5), c(0, 2),
    c(20, -40), c(35, -32), c(50, -25), c(55, -10), c(50, 0),
    c(53, 13), c(45, 10), c(37, 20), c(32, 30), c(16, 38), c(0, 38)
  )
)

check_area <- function(pop, map, continent) {
  if (is.null(pop)) return(NULL)

  # total population area
  pop_area <- area(pop)$area
  # population area overlapping a map
  map_area <- area(overlap(pop, map))
  # population area overlapping African continent
  continent_area <- area(overlap(pop, continent))

  # at least 50% of population's boundary be on land, and it must fall
  # on to the African continent itself
  if (continent_area == 0 || (map_area / pop_area) < 0.5)
    return(NULL)
  else
    return(pop)
}

filtered <- lapply(populations, check_area, map, continent) %>%
  Filter(Negate(is.null), .)

## ---- demes_africa------------------------------------------------------------
do.call(plot_map, filtered) + ggplot2::theme(legend.position = "none")

## ---- map_america-------------------------------------------------------------
xrange <- c(-90, -20)
yrange <- c(-58, 15)

map <- world(xrange = xrange, yrange = yrange, crs = "EPSG:31970")
plot_map(map)

## -----------------------------------------------------------------------------
# non-spatial ancestral population
p_anc <- population("p_anc", N = 1000, time = 1)

# spatial populations
p1 <- population("p1", N = 1000, time = 500, parent = p_anc, map = map, center = c(-75, 0), radius = 200e3)
p2 <- population("p2", N = 1000, time = 500, parent = p_anc, map = map, center = c(-60, 5), radius = 200e3)
p3 <- population("p3", N = 1000, time = 500, parent = p_anc, map = map, center = c(-65, -5), radius = 200e3)
p4 <- population("p4", N = 1000, time = 500, parent = p_anc, map = map, center = c(-60, -20), radius = 200e3)
p5 <- population("p5", N = 1000, time = 500, parent = p_anc, map = map, center = c(-65, -35), radius = 200e3)
p6 <- population("p6", N = 1000, time = 500, parent = p_anc, map = map, center = c(-69, -42), radius = 200e3)
p7 <- population("p7", N = 1000, time = 500, parent = p_anc, map = map, center = c(-51, -10), radius = 200e3)
p8 <- population("p8", N = 1000, time = 500, parent = p_anc, map = map, center = c(-45, -15), radius = 200e3)
p9 <- population("p9", N = 1000, time = 500, parent = p_anc, map = map, center = c(-71, -12), radius = 200e3)
p10 <- population("p10", N = 1000, time = 500, parent = p_anc, map = map, center = c(-55, -25), radius = 200e3)

## -----------------------------------------------------------------------------
gf <- list(
  gene_flow(p1, p2, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p2, p1, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p1, p3, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p3, p1, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p2, p3, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p3, p2, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p2, p7, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p7, p2, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p3, p7, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p7, p3, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p7, p8, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p8, p7, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p4, p7, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p7, p4, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p4, p5, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p5, p4, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p5, p6, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p6, p5, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p3, p4, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p4, p3, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p1, p9, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p9, p1, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p3, p9, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p9, p3, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p4, p9, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p9, p4, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p10, p4, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p4, p10, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p10, p8, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p8, p10, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p10, p5, 0.1, 1000, 2000, overlap = FALSE),
  gene_flow(p5, p10, 0.1, 1000, 2000, overlap = FALSE)
)

## -----------------------------------------------------------------------------
model <- compile_model(
  populations = list(p_anc, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10), gene_flow = gf,
  generation_time = 1, simulation_length = 5000,
  serialize = FALSE
)

## ---- model_america-----------------------------------------------------------
plot_model(model)

## ---- model_map_america-------------------------------------------------------
plot_map(model, gene_flow = TRUE)

## -----------------------------------------------------------------------------
ts <- msprime(model, sequence_length = 10e6, recombination_rate = 1e-8, random_seed = 42)

ts

