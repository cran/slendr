## ---- include = FALSE---------------------------------------------------------
env_present <- slendr:::is_slendr_env_present()
eval_chunk <- Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS"

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 80,
  eval = eval_chunk,
  fig.width = 8,
  fig.height = 6
)

set.seed(42)

## ---- message=FALSE-----------------------------------------------------------
library(ggplot2)
library(dplyr)

library(slendr)
init_env()

## ----introgression_model------------------------------------------------------
anc_all <- population("ancestor_all", time = 700e3, N = 10000, remove = 640e3)
afr <- population("AFR", parent = anc_all, time = 650e3, N = 10000)
anc_arch <- population("ancestor_archaics", parent = anc_all, time = 650e3, N = 10000, remove = 390e3)
nea <- population("NEA", parent = anc_arch, time = 400e3, N = 2000, remove = 30e3)
den <- population("DEN", parent = anc_arch, time = 400e3, N = 2000, remove = 30e3)
nonafr <- population("nonAFR", parent = afr, time = 100e3, N = 3000, remove = 39e3)
eur <- population("EUR", parent = nonafr, time = 45e3, N = 5000)
pap <- population("PAP", parent = nonafr, time = 45e3, N = 5000)

gf <- list(
  gene_flow(from = nea, to = nonafr, rate = 0.03, start = 55000, end = 50000),
  gene_flow(from = den, to = pap, rate = 0.07, start = 35000, end = 30000)
)

model <- compile_model(
  populations = list(anc_all, afr, anc_arch, nea, den, nonafr, eur, pap),
  gene_flow = gf,
  generation_time = 30,
  serialize = FALSE
)

plot_model(
  model, sizes = FALSE,
  order = c("AFR", "EUR", "nonAFR", "PAP", "ancestor_all", "DEN", "ancestor_archaics", "NEA")
)

## -----------------------------------------------------------------------------
samples <- schedule_sampling(model, times = 0, list(eur, 50), list(pap, 50))

ts <- msprime(model, sequence_length = 100e6, recombination_rate = 1e-8, samples = samples, random_seed = 42)

## -----------------------------------------------------------------------------
nea_tracts <- ts_tracts(ts, census = 55000)
den_tracts <- ts_tracts(ts, census = 35000)

tracts <- bind_rows(nea_tracts, den_tracts)

## -----------------------------------------------------------------------------
tracts

## -----------------------------------------------------------------------------
summary <- tracts %>%
  group_by(name, node_id, pop, source_pop) %>%
  summarise(prop = sum(length) / 100e6)

summary %>% group_by(pop, source_pop) %>% summarise(mean(prop)) %>% arrange(source_pop, pop)

## ---- anc_prop_summary--------------------------------------------------------
summary %>%
ggplot(aes(source_pop, prop, color = source_pop, fill = source_pop)) +
  geom_jitter() +
  coord_cartesian(ylim = c(0, 0.2)) +
  geom_hline(yintercept = c(0.03, 0.08), linetype = 2) +
  ylab("ancestry proportion") +
  facet_wrap(~ pop) +
  ggtitle("Ancestry proportions in each individual",
          "(vertical lines represent 3% and 7% baseline expectations")

## ----chrom_painting, fig.width=12, fig.height=10------------------------------
tracts %>%
mutate(chrom = paste(name, " (node", node_id, ")")) %>%
ggplot(aes(x = left, xend = right, y = chrom, yend = chrom, color = source_pop)) +
  geom_segment(linewidth = 3) +
  theme_minimal() +
  labs(x = "position [bp]", y = "haplotype") +
  ggtitle("True ancestry tracts along each chromosome") +
  theme(axis.text.y = element_blank(), panel.grid = element_blank()) +
  facet_grid(pop ~ ., scales = "free_y")

## -----------------------------------------------------------------------------
tracts %>%
  group_by(pop, source_pop) %>%
  summarise(mean(length))

## -----------------------------------------------------------------------------
m <- 0.03
t <- 52500 / 30
r <- 1e-8

mean_nea <- 1 / ((1 - m) * r * (t - 1))
mean_nea

## -----------------------------------------------------------------------------
m <- 0.07
t <- 37500 / 30
r <- 1e-8

mean_den <- 1 / ((1 - m) * r * (t - 1))
mean_den

## -----------------------------------------------------------------------------
expectation_df <- data.frame(
  pop = c("EUR", "PAP", "PAP"),
  source_pop = c("NEA", "NEA", "DEN"),
  length = c(mean_nea, mean_nea, mean_den)
)

## ----tract_lengths------------------------------------------------------------
p_densities <- tracts %>%
ggplot(aes(length, color = source_pop)) +
  geom_density() +
  geom_vline(data = expectation_df, aes(xintercept = length, color = source_pop),
             linetype = 2) +
  facet_wrap(~ pop) +
  ggtitle("Distribution of tract lengths per different ancestries")

cowplot::plot_grid(p_densities, p_densities + scale_x_log10(), nrow = 2)

## -----------------------------------------------------------------------------
sim_ts <- ts_load("/tmp/sim.trees")

squashed_tracts <- ts_tracts(sim_ts, census = 100.01, squashed = TRUE)

head(squashed_tracts)
tail(squashed_tracts)

## -----------------------------------------------------------------------------
full_tracts <- ts_tracts(sim_ts, census = 100.01, squashed = FALSE)

head(full_tracts)
tail(full_tracts)

