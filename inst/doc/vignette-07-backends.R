## ---- include = FALSE---------------------------------------------------------
env_present <- slendr:::is_slendr_env_present()

knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  dpi = 60,
  eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS"
)

## ---- message = FALSE---------------------------------------------------------
library(slendr)

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

init_env()

seed <- 42
set.seed(seed)

## -----------------------------------------------------------------------------
seq_len <- 100e6 # amount of sequence to simulate
rec_rate <- 1e-8 # uniform recombination rate
mut_rate <- 1e-8 # mutation rate

o <- population("outgroup", time = 1, N = 100)
c <- population("c", time = 2500, N = 100, parent = o)
a <- population("a", time = 3000, N = 100, parent = c)
b <- population("b", time = 3500, N = 100, parent = a)
x1 <- population("x1", time = 3800, N = 5000, parent = c)
x2 <- population("x2", time = 4000, N = 5000, parent = x1)

## -----------------------------------------------------------------------------
# no gene flow model
model_nogf <- compile_model(populations = list(a, b, x1, x2, c, o), generation_time = 1, simulation_length = 4500)

samples <- schedule_sampling(
  model_nogf, times = 4500,
  list(a, 1), list(b, 1), list(x1, 50), list(x2, 50), list(c, 1), list(o, 1)
)

# model with gene flow
gf <- gene_flow(from = b, to = x1, start = 4100, end = 4400, rate = 0.1)

model_gf <- compile_model(populations = list(a, b, x1, x2, c, o), gene_flow = gf, generation_time = 1, simulation_length = 4500)

samples <- schedule_sampling(
  model_gf, times = 4500,
  list(a, 1), list(b, 1), list(x1, 50), list(x2, 50), list(c, 1), list(o, 1)
)

## -----------------------------------------------------------------------------
plot_model(model_nogf, sizes = FALSE)
plot_model(model_gf, sizes = FALSE, proportions = TRUE)

## -----------------------------------------------------------------------------
# model without gene flow
slim_nogf <- slim(model_nogf, sequence_length = seq_len, recombination_rate = rec_rate, samples = samples, random_seed = seed)
msprime_nogf <- msprime(model_nogf, sequence_length = seq_len, recombination_rate = rec_rate, samples = samples, random_seed = seed)

# model with b -> x1 gene flow
slim_gf <- slim(model_gf, sequence_length = seq_len, recombination_rate = rec_rate, samples = samples, random_seed = seed)
msprime_gf <- msprime(model_gf, sequence_length = seq_len, recombination_rate = rec_rate, samples = samples, random_seed = seed)

## -----------------------------------------------------------------------------
# SLiM outputs -- we can use built-in slendr functions for those
slim_nogf <-
  slim_nogf %>%
  ts_recapitate(Ne = 10, recombination_rate = rec_rate, random_seed = seed) %>%
  ts_mutate(mut_rate, random_seed = seed)

slim_gf <-
  slim_gf %>%
  ts_recapitate(Ne = 10, recombination_rate = rec_rate, random_seed = seed) %>%
  ts_mutate(mut_rate, random_seed = seed)

# msprime outputs (note that recapitation and simplification doesn't make
# sense here because we already have fully coalesced genealogies for our
# individuals  of interest
msprime_nogf <- ts_mutate(msprime_nogf, mut_rate, random_seed = seed)
msprime_gf <- ts_mutate(msprime_gf, mut_rate, random_seed = seed)

## -----------------------------------------------------------------------------
slim_nogf
msprime_nogf

## -----------------------------------------------------------------------------
# extract vector of names of the "test individuals" in populations `x1` and `x2`
X <- ts_samples(slim_gf) %>% filter(pop %in% c("x1", "x2")) %>% pull(name)
X

## -----------------------------------------------------------------------------
# calculate f4-statistics on individuals of `x1` and `x2` populations using data
# from the two models (a model with no gene flow and a gene flow model) -- we use
# map_dfr to iterate across all individuals from `X_individuals` and binding all
# resulting data frames into a single data frame
df_slim_f4 <- rbind(
  map_dfr(X, ~ ts_f4(slim_nogf, "c_1", .x, "b_1", "outgroup_1")) %>% mutate(model = "no gene flow"),
  map_dfr(X, ~ ts_f4(slim_gf, "c_1", .x, "b_1", "outgroup_1")) %>% mutate(model = "gene flow")
) %>%
  select(X, f4, model) %>%
  mutate(simulator = "SLiM backend")

# compute the proportions of `b` ancestry in `x1` (expected 10%) and `x2`
# (expected 0% because this population did not receive any gene flow from `b`)
df_slim_f4ratio <- rbind(
  ts_f4ratio(slim_nogf, X, "a_1", "b_1", "c_1", "outgroup_1") %>% mutate(model = "no gene flow"),
  ts_f4ratio(slim_gf, X, "a_1", "b_1", "c_1", "outgroup_1") %>% mutate(model = "gene flow")
) %>%
  select(X, alpha, model) %>%
  mutate(simulator = "SLiM backend")

## -----------------------------------------------------------------------------
df_msprime_f4 <- rbind(
  map_dfr(X, ~ ts_f4(msprime_nogf, "c_1", .x, "b_1", "outgroup_1")) %>% mutate(model = "no gene flow"),
  map_dfr(X, ~ ts_f4(msprime_gf, "c_1", .x, "b_1", "outgroup_1")) %>% mutate(model = "gene flow")
) %>%
  select(X, f4, model) %>%
  mutate(simulator = "msprime backend")

# compute the proportions of `b` ancestry in `x1` (expected 10%) and `x2`
# (expected 0% because this population did not receive any gene flow from `b`)
df_msprime_f4ratio <- rbind(
  ts_f4ratio(msprime_nogf, X, "a_1", "b_1", "c_1", "outgroup_1") %>% mutate(model = "no gene flow"),
  ts_f4ratio(msprime_gf, X, "a_1", "b_1", "c_1", "outgroup_1") %>% mutate(model = "gene flow")
) %>%
  select(X, alpha, model) %>%
  mutate(simulator = "msprime backend")


## ---- msprime_slim_f4_distributions-------------------------------------------
df_f4 <- rbind(df_slim_f4, df_msprime_f4) %>%
  mutate(population = ifelse(grepl("x1_", X),
                             "x1 (received gene flow)",
                             "x2 (no gene flow)"))

ggplot(df_f4, aes(f4, fill = population)) +
  geom_histogram(bins = 50) +
  facet_grid(simulator ~ model) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(y = "number of individuals", x = "f4 statistic",
       title = "f4(c, x1 or x2; b, outgroup)",
       subtitle = "f4 ~0 is consistent with no gene flow, negative value indicates gene flow with 'b'") +
  theme(legend.position = "bottom")

## ---- msprime_slim_f4ratio_distributions--------------------------------------
df_f4ratio <- rbind(df_slim_f4ratio, df_msprime_f4ratio) %>%
  mutate(population = ifelse(grepl("x1_", X),
                             "x1 (received gene flow)",
                             "x2 (no gene flow)"))

ggplot(df_f4ratio, aes(alpha, fill = population)) +
  geom_histogram(bins = 30) +
  facet_grid(simulator ~ model) +
  geom_vline(xintercept = 0.1, linetype = 2) +
  labs(y = "number of individuals", x = "ancestry proportion (f4-ratio statistic)",
       title = "f4-ratio estimate of 'b' ancestry calculated from simulated data",
       subtitle = "f4-ratio = f4(a, outgroup; x1 or x2, c) / f4(a, outgroup; b, c)") +
  theme(legend.position = "bottom")

## -----------------------------------------------------------------------------
N <- 1000
N_factor <- 5 # by what factor should Ne change

seq_len <- 50e6
rec_rate <- 1e-8
mut_rate <- 1e-8

## -----------------------------------------------------------------------------
# constant Ne model
forward_const <- population("const", time = 1, N = N)

# decreasing step Ne model
forward_decr <- population("decr", time = 1, N = N, map = FALSE) %>%
  resize(time = 2000, N = N / N_factor, how = "step")

# increasing step Ne model
forward_incr <- population("inc", time = 1, N = N) %>%
  resize(time = 2000, N = N * N_factor, how = "step")

# exponential increase in size
forward_exp_incr <- population("exp_inc", time = 1, N = N) %>%
  resize(time = 2000, end = 3000, N = N * N_factor, how = "exponential")

# exponential decrease in size
forward_exp_decr <- population("exp_decr", time = 1, N = N) %>%
  resize(time = 2000, end = 3000, N = N / N_factor, how = "exponential")

## -----------------------------------------------------------------------------
# constant Ne model
backward_const <- population("const", time = 5000, N = N)

# decreasing step Ne model
backward_decr <- population("decr", time = 5000, N = N) %>%
  resize(time = 3000, N = N / N_factor, how = "step")

# increasing step Ne model
backward_incr <- population("inc", time = 5000, N = N) %>%
  resize(time = 3000, N = N * N_factor, how = "step")

# exponential increase in size
backward_exp_incr <- population("exp_inc", time = 5000, N = N) %>%
  resize(time = 3000, end = 2000, N = N * N_factor, how = "exponential")

# exponential decrease in size
backward_exp_decr <- population("exp_decr", time = 5000, N = N) %>%
  resize(time = 3000, end = 2000, N = N / N_factor, how = "exponential")

## -----------------------------------------------------------------------------
compile_run_afs <- function(model_name, pop, seed = 42) {
  # maximum length of the simulation (necessary for forward models which start
  # in generation 1)
  simulation_length <- 5000

  # define sampling times given the direction of time
  if (attr(pop, "history")[[1]]$time == 1) {
    sampling_time <- simulation_length
    direction <- "forward"
  } else {
    sampling_time <- 0
    direction <- "backward"
  }

  # compile model
  model <- compile_model(pop, generation_time = 15, direction = direction, simulation_length = simulation_length)

  samples <- schedule_sampling(model, times = sampling_time, list(pop, 50))

  # run the model in SLiM
  ts_slim <- slim(model, sequence_length = seq_len, recombination_rate = rec_rate,
                  samples = samples, random_seed = seed, verbose = FALSE)

  # run the same model in msprim
  ts_msprime <- msprime(model, sequence_length = seq_len, recombination_rate = rec_rate,
                        samples = samples, random_seed = seed, verbose = FALSE)

  # load the SLiM tree sequence
  ts_slim <- ts_recapitate(ts_slim, Ne = N, recombination_rate = rec_rate, random_seed = seed) %>%
    ts_mutate(mut_rate, random_seed = seed)

  # load the msprime tree sequence
  ts_msprime <- ts_mutate(ts_msprime, mut_rate, random_seed = seed)

  # compute the AFS from the SLiM and msprime tree sequences and bind the
  # results (derived allele counts per frequency bin) in a data frame
  msprime_afs <- ts_afs(ts_msprime, polarised = TRUE)[-1]
  slim_afs <- ts_afs(ts_slim, polarised = TRUE)[-1]

  rbind(
    data.frame(simulator = "msprime", model = model_name, f = msprime_afs),
    data.frame(simulator = "SLiM", model = model_name, f = slim_afs)
  ) %>%
    group_by(simulator, model) %>%
    mutate(n = 1:n(), direction = direction) %>%
    ungroup()
}

## -----------------------------------------------------------------------------
afs <- bind_rows(
  compile_run_afs("constant", forward_const),
  compile_run_afs("constant", backward_const),
  compile_run_afs("step contraction", forward_decr),
  compile_run_afs("step contraction", backward_decr),
  compile_run_afs("step increase", forward_incr),
  compile_run_afs("step increase", backward_incr),
  compile_run_afs("exponential decrease", forward_exp_decr),
  compile_run_afs("exponential decrease", backward_exp_decr),
  compile_run_afs("exponential increase", forward_exp_incr),
  compile_run_afs("exponential increase", backward_exp_incr)
) %>%
  mutate(model = factor(
    model, levels = c("step contraction", "constant", "step increase",
                      "exponential decrease", "exponential increase"))
  )

## ---- msprime_slim_afs--------------------------------------------------------
ggplot(afs, aes(n, f, color = direction, linetype = simulator)) +
  geom_line(stat = "identity") +
  facet_wrap(~ model) +
  labs(x = "number of derived alleles", y = "frequency",
       title = "Site frequency spectra obtained from five demographic models",
       subtitle = "Each model was specified in forward or backward direction of time and executed by
two different backend scripts in slendr (SLiM and msprime)") +
  guides(color = guide_legend("direction of\ntime in slendr"),
         linetype = guide_legend("slendr backend\nengine used")) +
  scale_linetype_manual(values = c(3, 2)) +
  scale_x_continuous(breaks = c(1, seq(20, 100, 20)), limits = c(1, 100)) +
  theme(legend.position = "bottom")

