## ---- include = FALSE---------------------------------------------------------
env_present <- slendr:::is_slendr_env_present()

knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  dpi = 80,
  eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") == ""
)

## -----------------------------------------------------------------------------
library(slendr)

# activate the internal Python environment needed for simulation and
# tree-sequence processing
init_env()

## ---- world_zoom, results = FALSE---------------------------------------------
map <- world(
  xrange = c(-13, 70), # min-max longitude
  yrange = c(18, 65),  # min-max latitude
  crs = "EPSG:3035"    # coordinate reference system (CRS) for West Eurasia
)

## -----------------------------------------------------------------------------
map

## -----------------------------------------------------------------------------
africa <- region(
  "Africa", map,
  polygon = list(c(-18, 20), c(38, 20), c(30, 33),
                 c(20, 33), c(10, 38), c(-6, 35))
)
europe <- region(
  "Europe", map,
  polygon = list(
    c(-8, 35), c(-5, 36), c(10, 38), c(20, 35), c(25, 35),
    c(33, 45), c(20, 58), c(-5, 60), c(-15, 50)
  )
)
anatolia <- region(
  "Anatolia", map,
  polygon = list(c(28, 35), c(40, 35), c(42, 40),
                 c(30, 43), c(27, 40), c(25, 38))
)

## -----------------------------------------------------------------------------
anatolia

## -----------------------------------------------------------------------------
class(anatolia)

## ----plot_europe_anatolia-----------------------------------------------------
plot_map(africa, europe, anatolia, title = "Geographic regions")

## -----------------------------------------------------------------------------
all(attr(europe, "map") == map)
all(attr(anatolia, "map") == map)

## ----plot_afr-----------------------------------------------------------------
afr <- population("AFR", time = 52000, N = 3000, map = map, polygon = africa)

plot_map(afr)

## ---- message = FALSE---------------------------------------------------------
ooa <- population(
  "OOA", parent = afr, time = 51000, N = 500, remove = 25000,
  center = c(33, 30), radius = 400e3
)

## ----plot_ooa-----------------------------------------------------------------
plot_map(ooa, intersect = TRUE, title = "'Intersected' population range")

## ---- message = FALSE---------------------------------------------------------
ooa <- ooa %>% move(
  trajectory = list(c(40, 30), c(50, 30), c(60, 40)),
  start = 50000, end = 40000
)

## -----------------------------------------------------------------------------
ooa

## ----plot_ooa_migration-------------------------------------------------------
plot_map(ooa, title = "Intermediate migration maps")

## ----plot_ehg-----------------------------------------------------------------
ehg <- population(
  "EHG", parent = ooa, time = 28000, N = 1000, remove = 6000,
  polygon = list(
    c(26, 55), c(38, 53), c(48, 53), c(60, 53),
    c(60, 60), c(48, 63), c(38, 63), c(26, 60))
)

## ----plot_eur-----------------------------------------------------------------
eur <- population( # European population
  name = "EUR", parent = ehg, time = 25000, N = 2000,
  polygon = europe
)

## ---- message = FALSE---------------------------------------------------------
ana <- population( # Anatolian farmers
  name = "ANA", time = 28000, N = 3000, parent = ooa, remove = 4000,
  center = c(34, 38), radius = 500e3, polygon = anatolia
) %>%
  expand_range( # expand the range by 2.500 km
    by = 2500e3, start = 10000, end = 7000,
    polygon = join(europe, anatolia)
  )

## -----------------------------------------------------------------------------
ana

## ----plot_ana-----------------------------------------------------------------
plot_map(ana, title = "Anatolian expansion into Europe")

## ----plot_ana_raw, eval = FALSE-----------------------------------------------
#  plot_map(ana, title = "Anatolian expansion into Europe (not intersected)", intersect = FALSE)

## ----plot_yam_migr------------------------------------------------------------
yam <- population( # Yamnaya steppe population
  name = "YAM", time = 7000, N = 500, parent = ehg, remove = 2500,
  polygon = list(c(26, 50), c(38, 49), c(48, 50),
                 c(48, 56), c(38, 59), c(26, 56))
) %>%
  move(
    trajectory = c(15, 50),
    start = 5000, end = 3000, snapshots = 8
  )

plot_map(yam)

## ----plot_maps----------------------------------------------------------------
plot_map(afr, ooa, ehg, eur, ana, yam)

## ----eval = FALSE-------------------------------------------------------------
#  gf <- gene_flow(from = eur, to = afr, rate = 0.1, start = 20000, end = 15000)

## -----------------------------------------------------------------------------
gf <- list(
  gene_flow(from = ana, to = yam, rate = 0.5, start = 6500, end = 6400, overlap = FALSE),
  gene_flow(from = ana, to = eur, rate = 0.5, start = 8000, end = 6000),
  gene_flow(from = yam, to = eur, rate = 0.75, start = 4000, end = 3000)
)

## -----------------------------------------------------------------------------
gf

## -----------------------------------------------------------------------------
model_dir <- paste0(tempfile(), "_tutorial-model")

model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam), # populations defined above
  gene_flow = gf, # gene-flow events defined above
  generation_time = 30,
  resolution = 10e3, # resolution in meters per pixel
  competition = 130e3, mating = 100e3, # spatial interaction in SLiM
  dispersal = 70e3, # how far will offspring end up from their parents
  path = model_dir
)

## -----------------------------------------------------------------------------
list.files(model_dir, pattern = "*.jpg")

## -----------------------------------------------------------------------------
read.table(file.path(model_dir, "populations.tsv"), header = TRUE)

## -----------------------------------------------------------------------------
read.table(file.path(model_dir, "geneflow.tsv"), header = TRUE)

## -----------------------------------------------------------------------------
loaded_model <- read_model(model_dir)

## ----plot_model, fig.width = 8, fig.height = 7--------------------------------
plot_model(model, proportions = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  explore_model(model)

## -----------------------------------------------------------------------------
ts <- slim(model, sequence_length = 100000, recombination_rate = 1e-8)
ts

