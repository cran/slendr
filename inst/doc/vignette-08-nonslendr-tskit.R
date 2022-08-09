## ---- include = FALSE---------------------------------------------------------
env_present <- slendr:::check_env_present()

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  dpi = 80,
  eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS"
)

nonspatial_slim <- tempfile()
nonspatial_trees_file <- tempfile()

spatial_slim <- tempfile()
spatial_trees_file <- tempfile()

msprime_trees_file <- tempfile()

Sys.setenv(NONSPATIAL_SLIM = nonspatial_slim)
Sys.setenv(NONSPATIAL_TREES = nonspatial_trees_file)

Sys.setenv(SPATIAL_SLIM = spatial_slim)
Sys.setenv(SPATIAL_TREES = spatial_trees_file)

Sys.setenv(MSPRIME_TREES = msprime_trees_file)

## ---- message = FALSE---------------------------------------------------------
library(slendr)

library(dplyr)
library(magrittr)
library(ggplot2)

SEED <- 42
set.seed(SEED)

## -----------------------------------------------------------------------------
ts <- ts_load(nonspatial_trees_file) %>%
  ts_simplify() %>%
  ts_mutate(mutation_rate = 1e-7, random_seed = SEED)

## -----------------------------------------------------------------------------
data <- ts_nodes(ts) %>% dplyr::filter(sampled)
data

## -----------------------------------------------------------------------------
sample_sets <- split(data$node_id, data$pop)

# compute nucleotide diversity in each population
# (any other ts_*() tskit R interface function should work)
ts_diversity(ts, sample_sets)

## -----------------------------------------------------------------------------
samples <- sample(data$node_id, 10)
ts_small <- ts_simplify(ts, simplify_to = samples)

# extract the 42nd tree in the genealogy to an R 'phylo' format
tree <- ts_phylo(ts_small, 42)
tree

## ---- nonslendr_tree, eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS" && Sys.getenv("R_HAS_GGTREE") == TRUE----
#  library(ggtree)
#  
#  labels <- ts_nodes(tree) %>% select(node = phylo_id, tskit_id = node_id)
#  
#  ggtree(tree, branch.length = "none") %<+% labels +
#    geom_label(aes(label = tskit_id))

## ---- eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS" && Sys.getenv("R_HAS_GGTREE") != TRUE----
library(ape)
plot(tree, show.tip.label = FALSE)
nodelabels()
tiplabels()

## ---- include = FALSE---------------------------------------------------------
reticulate::py_run_string(sprintf("import msprime; msprime.simulate(100).dump('%s')", msprime_trees_file))

## -----------------------------------------------------------------------------
ts <- ts_load(msprime_trees_file)

ts_nodes(ts)

## -----------------------------------------------------------------------------
ts <- ts_load(spatial_trees_file) %>% ts_simplify()

## -----------------------------------------------------------------------------
data <- ts_nodes(ts)
data

## ---- nonslendr_locations-----------------------------------------------------
ggplot() + geom_sf(data = data, aes(color = time), alpha = 0.5)

## ---- nonslendr_ancestries----------------------------------------------------
ancestral_links <- ts_ancestors(ts, 0)

ggplot() +
    geom_sf(data = ancestral_links, size = 0.5, aes(alpha = parent_time)) +
    geom_sf(data = sf::st_set_geometry(ancestral_links, "parent_location"), aes(color = parent_time)) +
    geom_sf(data = data[data$node_id == 0, ], size = 3, color = "red")

