## ----include = FALSE----------------------------------------------------------
env_present <- slendr:::is_slendr_env_present()

knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  eval = env_present
)

## -----------------------------------------------------------------------------
library(slendr)

init_env()

# African ancestral population
afr <- population("AFR", time = 52000, N = 3000)

# first migrants out of Africa
ooa <- population("OOA", parent = afr, time = 51000, N = 500, remove = 25000)

# Eastern hunter-gatherers
ehg <- population("EHG", parent = ooa, time = 28000, N = 1000, remove = 6000)

# European population
eur <- population("EUR", parent = ehg, time = 25000, N = 2000) %>%
  resize(N = 10000, how = "exponential", time = 5000, end = 0)

# Anatolian farmers
ana <- population("ANA", time = 28000, N = 3000, parent = ooa, remove = 4000)

# Yamnaya steppe population
yam <- population("YAM", time = 7000, N = 500, parent = ehg, remove = 2500)

## -----------------------------------------------------------------------------
gf <- list(
  gene_flow(from = ana, to = yam, rate = 0.5, start = 6500, end = 6400),
  gene_flow(from = ana, to = eur, rate = 0.5, start = 8000, end = 6000),
  gene_flow(from = yam, to = eur, rate = 0.75, start = 4000, end = 3000)
)

## -----------------------------------------------------------------------------
model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30
)

## ----non-spatial_graph, fig.width = 6, fig.height = 6, dpi = 60---------------
plot_model(model)

## ----eval = FALSE-------------------------------------------------------------
#  ts_slim <- slim(model, sequence_length = 100000, recombination_rate = 0)

## ----eval = FALSE-------------------------------------------------------------
#  ts_msprime <- msprime(model, sequence_length = 100000, recombination_rate = 0)

