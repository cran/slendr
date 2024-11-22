## ----include = FALSE----------------------------------------------------------
env_present <- slendr:::is_slendr_env_present()

#if (Sys.getenv("RUNNER_OS") == "Windows") {
#  bash_path <- "D:/a/_temp/msys64/usr/bin/bash.exe"
#} else {
#  bash_path <- NULL
#}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  dpi = 60,
  eval = slendr:::is_slim_present() && env_present && Sys.getenv("RUNNER_OS") != "Windows"#,
  #engine.path = list(bash = bash_path)
)

nonspatial_slim <- normalizePath(tempfile(), winslash = "/", mustWork = FALSE)
nonspatial_trees_file <-normalizePath(tempfile(), winslash = "/", mustWork = FALSE)

spatial_slim <- normalizePath(tempfile(), winslash = "/", mustWork = FALSE)
spatial_trees_file <- normalizePath(tempfile(), winslash = "/", mustWork = FALSE)

msprime_trees_file <- normalizePath(tempfile(), winslash = "/", mustWork = FALSE)

Sys.setenv(NONSPATIAL_SLIM = nonspatial_slim)
Sys.setenv(NONSPATIAL_TREES = nonspatial_trees_file)

Sys.setenv(SPATIAL_SLIM = spatial_slim)
Sys.setenv(SPATIAL_TREES = spatial_trees_file)

Sys.setenv(MSPRIME_TREES = msprime_trees_file)

## ----message = FALSE----------------------------------------------------------
library(slendr)

library(dplyr)
library(magrittr)
library(ggplot2)

init_env()

SEED <- 42
set.seed(SEED)

## -----------------------------------------------------------------------------
ts <- ts_read(nonspatial_trees_file) %>%
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
tree <- ts_phylo(ts_small, 42 - 1)
tree

## ----nonslendr_tree, eval = slendr:::is_slim_present() && env_present && Sys.getenv("R_HAS_GGTREE") == TRUE && Sys.getenv("RUNNER_OS") != "Windows"----
#  library(ggtree)
#  
#  labels <- ts_nodes(tree) %>% select(node = phylo_id, tskit_id = node_id)
#  
#  ggtree(tree, branch.length = "none") %<+% labels +
#    geom_label(aes(label = tskit_id))

## ----eval = slendr:::is_slim_present() && env_present && Sys.getenv("R_HAS_GGTREE") != TRUE && Sys.getenv("RUNNER_OS") != "Windows"----
library(ape)
plot(tree, show.tip.label = FALSE)
nodelabels()
tiplabels()

## ----include = FALSE----------------------------------------------------------
reticulate::py_run_string(sprintf("import msprime; msprime.simulate(100).dump('%s')", msprime_trees_file))

## -----------------------------------------------------------------------------
ts <- ts_read(msprime_trees_file)

ts_nodes(ts)

## -----------------------------------------------------------------------------
ts <- ts_read(spatial_trees_file) %>% ts_simplify()

## -----------------------------------------------------------------------------
data <- ts_nodes(ts)
data

## ----nonslendr_locations------------------------------------------------------
ggplot() + geom_sf(data = data, aes(color = time), alpha = 0.5)

## ----nonslendr_ancestries-----------------------------------------------------
ancestral_links <- ts_ancestors(ts, 0)

ggplot() +
    geom_sf(data = ancestral_links, size = 0.5, aes(alpha = parent_time)) +
    geom_sf(data = sf::st_set_geometry(ancestral_links, "parent_location"), aes(color = parent_time)) +
    geom_sf(data = data[data$node_id == 0, ], size = 3, color = "red")

