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

set.seed(314159)

## ----abstract_world-----------------------------------------------------------
map <- world(xrange = c(0, 3500), yrange = c(0, 700), landscape = "blank")

## ----pop_ranges---------------------------------------------------------------
N <- 3000; y <- 350; r = 240

p0 <- population("pop0", time = 1, N = N, map = map, center = c(250, y), radius = r)
p1 <- population("pop1", time = 1, N = N, map = map, center = c(750, y), radius = r)
p2 <- population("pop2", time = 1, N = N, map = map, center = c(1250, y), radius = r)
p3 <- population("pop3", time = 1, N = N, map = map, center = c(1750, y), radius = r)
p4 <- population("pop4", time = 1, N = N, map = map, center = c(2250, y), radius = r)
p5 <- population("pop5", time = 1, N = N, map = map, center = c(2750, y), radius = r)
p6 <- population("pop6", time = 1, N = N, map = map, center = c(3250, y), radius = r)

plot_map(p0, p1, p2, p3, p4, p5, p6)

## -----------------------------------------------------------------------------
p1 <- set_dispersal(p1, time = 100, competition = 80)
p2 <- set_dispersal(p2, time = 200, competition = 130)
p3 <- set_dispersal(p3, time = 300, competition = 170)
p4 <- set_dispersal(p4, time = 400, competition = 220)
p5 <- set_dispersal(p5, time = 500, competition = 300)
p6 <- set_dispersal(p6, time = 600, competition = 380)

## -----------------------------------------------------------------------------
model <- compile_model(
  populations = list(p0, p1, p2, p3, p4, p5, p6),
  generation_time = 1, resolution = 1, simulation_length = 1000,
  competition = 10, mating = 10, dispersal = 5,
  path = paste0(tempfile(), "_spatial-interactions")
)

## -----------------------------------------------------------------------------
locations_file <- tempfile(fileext = ".gz")
ts <- slim(model, sequence_length = 1, recombination_rate = 0, locations = locations_file,
           burnin = 100, method = "batch", verbose = FALSE, random_seed = 314159)

# get a summary of the simulated tree-sequence object
ts

## ----plot_gif_interactions, message = FALSE, eval = FALSE---------------------
#  animate_model(model, locations_file, steps = 80, width = 500, height = 200)

## ----plot_interactions, eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS"----
library(ggplot2)

locations <- ts_nodes(ts) %>% dplyr::filter()
ggplot() + geom_sf(data = locations, aes(color = pop), size = 0.75, alpha = 0.5) + theme_minimal()

