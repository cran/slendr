## ---- include = FALSE---------------------------------------------------------
env_present <- slendr:::is_slendr_env_present()

knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  dpi = 80,
  eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS" && Sys.which("qpDstat") != ""
)

## -----------------------------------------------------------------------------
library(slendr)

## ---- eval = FALSE------------------------------------------------------------
#  setup_env()

## -----------------------------------------------------------------------------
init_env()

## -----------------------------------------------------------------------------
check_env()

## -----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)

set.seed(314159)

# create the ancestor of everyone and a chimpanzee outgroup
# (we set both N = 1 to reduce the computational time for this model)
chimp <- population("CH", time = 6.5e6, N = 1000)

# two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 10000)
eur <- population("EUR", parent = afr, time = 70e3, N = 5000)

# Neanderthal population splitting at 600 ky ago from modern humans
# (becomes extinct by 40 ky ago)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

# 3% Neanderthal introgression into Europeans between 55-50 ky ago
gf <- gene_flow(from = nea, to = eur, rate = 0.03, start = 55000, end = 45000)

model <- compile_model(
  populations = list(chimp, nea, afr, eur), gene_flow = gf,
  generation_time = 30,
  path = paste0(tempfile(), "_introgression")
)

## ---- introgression_graph, fig.width = 8, fig.height = 4----------------------
cowplot::plot_grid(
  plot_model(model, sizes = FALSE),
  plot_model(model, sizes = FALSE, log = TRUE),
  nrow = 1
)

## -----------------------------------------------------------------------------
nea_samples <- schedule_sampling(model, times = c(70000, 40000), list(nea, 1))
nea_samples

## -----------------------------------------------------------------------------
present_samples <- schedule_sampling(model, times = 0, list(chimp, 1), list(afr, 5), list(eur, 10))
present_samples

## -----------------------------------------------------------------------------
emh_samples <- schedule_sampling(model, times = runif(n = 40, min = 10000, max = 40000), list(eur, 1))
emh_samples

## ---- eval = FALSE------------------------------------------------------------
#  # this attempts to sample a Neanderthal individual at a point when Neanderthals
#  # are already extinct, resulting in an error
#  schedule_sampling(model, times = 10000, list(nea, 1), strict = TRUE)

## -----------------------------------------------------------------------------
ts <- msprime(
  model, sequence_length = 100e6, recombination_rate = 1e-8,
  samples = rbind(nea_samples, present_samples, emh_samples),
  random_seed = 314159, verbose = TRUE
)

ts

## -----------------------------------------------------------------------------
output_file <- tempfile()

ts <- msprime(
  model, sequence_length = 100e6, recombination_rate = 1e-8,
  samples = rbind(nea_samples, present_samples, emh_samples),
  output = output_file, random_seed = 314159
)

output_file

## -----------------------------------------------------------------------------
ts <- ts_load(output_file, model)
ts

## -----------------------------------------------------------------------------
ts_simplify(ts)

## -----------------------------------------------------------------------------
ts_small <- ts_simplify(ts, simplify_to = c("CH_1", "NEA_1", "NEA_2", "AFR_1", "AFR_2", "EUR_20", "EUR_50"))
ts_small

## ---- eval = FALSE------------------------------------------------------------
#  ts <- ts_recapitate(ts, recombination_rate = 1e-8, Ne = 10000)

## -----------------------------------------------------------------------------
ts_coalesced(ts)

## -----------------------------------------------------------------------------
ts <- ts_mutate(ts, mutation_rate = 1e-8, random_seed = 314159)
ts

## -----------------------------------------------------------------------------
# extract the 42nd tree in the tree sequence
tree <- ts_phylo(ts_small, 42)

## -----------------------------------------------------------------------------
tree

## ---- plot_ggtree, eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS" && Sys.getenv("R_HAS_GGTREE") == TRUE----
#  library(ggtree)
#  
#  ggtree(tree) +
#    geom_point2(aes(subset = !isTip)) + # points for internal nodes
#    geom_tiplab() + # sample labels for tips
#    hexpand(0.1)    # make more space for the tip labels

## ---- plot_apetree, eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS" && Sys.getenv("R_HAS_GGTREE") != TRUE----
library(ape)
plot(tree)
nodelabels()

## ---- eval = FALSE------------------------------------------------------------
#  # iterate over all trees in the tree-sequence and check if each
#  # has only one root (i.e. is fully coalesced) - note that Python
#  # lists are 0-based, which is something we need to take care of
#  all(sapply(seq_len(ts$num_trees)[1:100],
#             function(i) ts$at_index(i - 1)$num_roots == 1))

## -----------------------------------------------------------------------------
names(ts)

## -----------------------------------------------------------------------------
# f2 is a measure of the branch length connecting A and B
ts_f2(ts, A = "EUR_1", B = "AFR_1")

# f4 is a measure of the drift shared between A and B after their split from C
ts_f3(ts, A = "EUR_1", B = "AFR_1", C = "CH_1")

# this value should be very close to zero (no introgression in Africans)
ts_f4(ts, "AFR_1", "AFR_2", "NEA_1", "CH_1", mode = "branch")

# this value should be significantly negative (many more ABBA sites
# compared to BABA site due to the introgression into Europeans)
ts_f4(ts, "AFR_1", "EUR_1", "NEA_1", "CH_1", mode = "branch")

## -----------------------------------------------------------------------------
ts_samples(ts)

## -----------------------------------------------------------------------------
# first get a table of simulated African and European individuals in the tree-sequence
inds <- ts_samples(ts) %>% dplyr::filter(pop %in% c("AFR", "EUR"))

# estimate the amounts of Neanderthal ancestry in these individuals and add
# these values to the table
inds$ancestry <- ts_f4ratio(ts, X = inds$name, "NEA_1", "NEA_2", "AFR_1", "CH_1")$alpha

## ---- ts_nea_distributions----------------------------------------------------
ggplot(inds, aes(pop, ancestry, fill = pop)) +
  geom_boxplot() +
  geom_jitter() +
  labs(y = "Neanderthal ancestry proportion", x = "") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.1))

## ---- ts_nea_trajectory-------------------------------------------------------
dplyr::filter(inds, pop == "EUR") %>%
  ggplot(aes(time, ancestry)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = 2, color = "red", linewidth = 0.5) +
  xlim(40000, 0) + coord_cartesian(ylim = c(0, 0.1)) +
  labs(x = "time [years ago]", y = "Neanderthal ancestry proportion")

## -----------------------------------------------------------------------------
snps <- ts_eigenstrat(ts, prefix = file.path(tempdir(), "eigenstrat", "data"))

## ---- eval = Sys.which("qpDstat") != "" && env_present && Sys.getenv("R_HAS_GGTREE") == TRUE----
#  library(admixr)
#  
#  f4ratio(data = snps, X = c("EUR_1", "EUR_2", "AFR_2"),
#          A = "NEA_1", B = "NEA_2", C = "AFR_1", O = "CH_1")

## ---- ts_vs_admixr, eval = Sys.which("qpDstat") != "" && env_present && Sys.getenv("R_HAS_GGTREE") == TRUE----
#  europeans <- inds[inds$pop == "EUR", ]$name
#  
#  # tskit result
#  result_ts <- ts_f4ratio(ts, X = europeans, A = "NEA_1", B = "NEA_2", C = "AFR_1", O = "CH_1") %>% select(alpha_ts = alpha)
#  
#  # result obtained by admixr/ADMIXTOOLS
#  result_admixr <- f4ratio(snps, X = europeans, A = "NEA_1", B = "NEA_2", C = "AFR_1", O = "CH_1") %>% select(alpha_admixr = alpha)
#  
#  bind_cols(result_admixr, result_ts) %>%
#    ggplot(aes(alpha_ts, alpha_admixr)) +
#    geom_point() +
#    geom_abline(slope = 1, linetype = 2, color = "red", linewidth = 0.5) +
#    labs(x = "f4-ratio statistic calculated with admixr/ADMIXTOOLS",
#         y = "f4-ratio statistic calculated with tskit")

## -----------------------------------------------------------------------------
ts_vcf(ts, path = file.path(tempdir(), "output.vcf.gz"))

## -----------------------------------------------------------------------------
ts_vcf(ts, path = file.path(tempdir(), "output_subset.vcf.gz"),
       individuals = c("CH_1", "NEA_1", "EUR_1", "AFR_1"))

## -----------------------------------------------------------------------------
ts_fst(ts, sample_sets = list(afr = c("AFR_1", "AFR_2", "AFR_3"), eur = c("EUR_1", "EUR_2")))

## -----------------------------------------------------------------------------
ts_fst(ts, sample_sets = list(c("AFR_1", "AFR_2", "AFR_3"), c("EUR_1", "EUR_2")))

## -----------------------------------------------------------------------------
ts_fst(ts, sample_sets = list(afr = c("AFR_1", "AFR_2", "AFR_3"),
                              eur = c("EUR_1", "EUR_2"),
                              nea = c("NEA_1", "NEA_2")))

## -----------------------------------------------------------------------------
# define breakpoints between 20 windows
breakpoints <- seq(0, ts$sequence_length, length.out = 21)

# calculate window-based Fst statistic
win_fst <- ts_fst(
  ts, windows = breakpoints,
  sample_sets = list(afr = c("AFR_1", "AFR_2", "AFR_3"),
                     eur = c("EUR_1", "EUR_2"),
                     nea = c("NEA_1", "NEA_2"))
)

# we get 20 values for each parwise calculation
win_fst

## -----------------------------------------------------------------------------
win_fst[1, ]$Fst

## -----------------------------------------------------------------------------
ts_tajima(ts, list(afr = c("AFR_1", "AFR_2", "AFR_3"), eur = c("EUR_1", "EUR_2")))

## -----------------------------------------------------------------------------
ts_tajima(ts, list(afr = c("AFR_1", "AFR_2"), eur = c("EUR_1", "EUR_2")), windows = breakpoints)

## -----------------------------------------------------------------------------
# get sampled individuals from all populations
sample_sets <- ts_samples(ts) %>%
  split(., .$pop) %>%
  lapply(function(pop) pop$name)

sample_sets

## -----------------------------------------------------------------------------
ts_diversity(ts, sample_sets) %>% dplyr::arrange(diversity)

## -----------------------------------------------------------------------------
ts_divergence(ts, sample_sets) %>% arrange(divergence)

