## ----include = FALSE----------------------------------------------------------
run_vignette <- Sys.getenv("RUNNER_OS") == "" && slendr::check_dependencies(python = TRUE, slim = TRUE, quit = FALSE)

knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  dpi = 60,
  eval = run_vignette
)

RERUN <- FALSE

## ----collapse = TRUE, message = FALSE-----------------------------------------
library(slendr)
init_env()

library(dplyr)
library(ggplot2)
library(readr)

seed <- 42
set.seed(seed)

## ----eval=FALSE---------------------------------------------------------------
# afr <- population("AFR", time = 100000, N = 20000)
# eur <- population("EUR", time = 60000, N = 2000, parent = afr)
# 
# # <... compile model and simulate a tree sequence `ts` ...>

## ----eval=FALSE---------------------------------------------------------------
# # compute heterozygosity in the individual "EUR"
# ts_diversity(ts, "EUR_1")
# 
# # compute genetic divergence between selected Africans and Europeans
# afr_samples <- c("AFR_1", "AFR_2", "AFR_3")
# eur_samples <- c("EUR_1", "EUR_2", "EUR_3")
# ts_divergence(ts, list(afr = afr_samples, eur = eur_samples))

## ----eval=FALSE---------------------------------------------------------------
# afr <- population("AFR", time = 90000, N = 20000)
# eur <- population("EUR", time = 60000, N = 2000, parent = afr)

## ----eval=FALSE---------------------------------------------------------------
# model <- compile_model(populations = list(afr, eur), generation_time = 30)
# 
# slim(model, sequence_length = 100000, recombination_rate = 1e-8, method = "gui")

## ----eval=FALSE---------------------------------------------------------------
# pop <- population("pop", time = 1, N = 1000)
# simple_model <- compile_model(pop, generation_time = 1, simulation_length = 1000)
# 
# slim(simple_model, sequence_length = 1000, recombination_rate = 0, method = "gui")

## ----eval=FALSE---------------------------------------------------------------
# pop <- population("pop", time = 1, N = 1000)
# simple_model <- compile_model(pop, generation_time = 1, simulation_length = 1000)
# 
# slim(simple_model, sequence_length = 1000, recombination_rate = 0, burnin = 100, method = "gui")

## -----------------------------------------------------------------------------
library(slendr)
init_env()

# African ancestral population
afr <- population("AFR", time = 90000, N = 3000)

# first migrants out of Africa
ooa <- population("OOA", parent = afr, time = 60000, N = 500, remove = 23000) %>%
  resize(N = 2000, time = 40000, how = "step")

# Eastern hunter-gatherers
ehg <- population("EHG", parent = ooa, time = 28000, N = 1000, remove = 6000)

# European population
eur <- population("EUR", parent = ehg, time = 25000, N = 5000)

# Anatolian farmers
ana <- population("ANA", time = 28000, N = 3000, parent = ooa, remove = 4000)

# Yamnaya steppe population
yam <- population("YAM", time = 8000, N = 500, parent = ehg, remove = 2500)

# define gene-flow events
gf <- list(
  gene_flow(from = ana, to = yam, rate = 0.4, start = 7900, end = 7800),
  gene_flow(from = ana, to = eur, rate = 0.5, start = 6000, end = 5000),
  gene_flow(from = yam, to = eur, rate = 0.65, start = 4000, end = 3500)
)

## ----model_plot, echo=FALSE---------------------------------------------------
model_for_plot <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30
)

plot_model(model_for_plot, order = c("EUR", "EHG", "YAM", "ANA", "OOA", "AFR"), proportions = TRUE)

## -----------------------------------------------------------------------------
extension_path <- system.file("extdata", "extension_trajectory.txt", package = "slendr")

## ----echo=FALSE---------------------------------------------------------------
Sys.setenv(EXTENSION1 = extension_path)

## -----------------------------------------------------------------------------
model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30,
  extension = extension_path  # <--- include the SLiM extension snippet
)

## ----eval=FALSE---------------------------------------------------------------
# slim(model, sequence_length = 1e6, recombination_rate = 0, method = "gui")

## ----echo=FALSE---------------------------------------------------------------
slim(model, sequence_length = 1e6, recombination_rate = 0, path = tempdir())

## -----------------------------------------------------------------------------
extension_path <- system.file("extdata", "extension_trajectory_params.txt", package = "slendr")

## ----echo=FALSE---------------------------------------------------------------
Sys.setenv(EXTENSION = extension_path)

## -----------------------------------------------------------------------------
extension <- substitute_values(
  extension_path,
  s = 0.1, onset_time = 15000,
  origin_pop = "EUR", target_pop = "EUR"
)

## ----eval = FALSE-------------------------------------------------------------
# extension <- substitute_values(
#   extension_path,
#   onset_time = 15000,
#   origin_pop = "EUR", target_pop = "EUR"
# )
# 
# # Error: The extension script contains the following unsubstituted patterns: {{s}}

## -----------------------------------------------------------------------------
model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf, generation_time = 30,
  extension = extension
)

## ----echo = FALSE-------------------------------------------------------------
slim(model, sequence_length = 1e6, recombination_rate = 1e-8, path = tempdir())

## ----eval = FALSE-------------------------------------------------------------
# slim(model, sequence_length = 1e6, recombination_rate = 0, path = tempdir(), random_seed = 42)

## -----------------------------------------------------------------------------
run_model <- function(origin_pop, onset_time) {
  extension <- substitute_values(
    extension_path,
    s = 0.1, onset_time = onset_time,
    origin_pop = origin_pop, target_pop = "EUR"
  )

  model <- compile_model(
    populations = list(afr, ooa, ehg, eur, ana, yam),
    gene_flow = gf, generation_time = 30,
    extension = extension
  )

  slim(model, sequence_length = 1e6, recombination_rate = 0,
       path = tempdir(), random_seed = 42)
}

run_model(origin_pop = "EUR", onset_time = 15000)
run_model(origin_pop = "ANA", onset_time = 15000)
run_model(origin_pop = "EHG", onset_time = 15000)
run_model(origin_pop = "YAM", onset_time = 8000)

## ----trajectories-------------------------------------------------------------
load_traj <- function(origin_pop) {
  df <- read.table(paste0(tempdir(), "/traj_EUR_", origin_pop, ".tsv"), header = TRUE)
  df$origin <- origin_pop
  df$target <- "EUR"
  df
}

traj <- rbind(load_traj("EUR"), load_traj("ANA"), load_traj("EHG"), load_traj("YAM"))

library(ggplot2)

ggplot(traj) +
  geom_line(aes(time, freq_target, linetype = "EUR"), color = "black") +
  geom_line(aes(time, freq_origin, color = origin), linetype = "dashed") +
  xlim(15000, 0) +
  labs(title = "Allele frequency in EUR given the origin in another population",
       x = "years before present", y = "allele frequency",
       color = "frequency\nin original\npopulation",
       linetype = "frequency\nin target\npopulation") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  facet_wrap(~ origin)

## -----------------------------------------------------------------------------
# extension "template" provided as a single string (this contains the same code
# as the script used just above, except specified directly in R)
extension_template <- r"(
// Because we want to simulate non-neutral evolution, we have to provide a
// custom initialization callback -- slendr will use it to replace its default
// neutral genomic architecture (i.e. the initialize() {...} callback it uses
// by default for neutral simulations). Note that we can refer to slendr's
// constants SEQUENCE_LENGTH and RECOMBINATION_RATE, which will carry values
// passed through from R via slendr's slim() R function.
initialize() {
    initializeMutationType("m1", 0.5, "f", 0.0);

    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, SEQUENCE_LENGTH - 1);

    initializeMutationRate(0);
    initializeRecombinationRate(RECOMBINATION_RATE);
}

// Define model constants (to be substituted) all in one place
// (each {{placeholder}} will be replaced by a value passed from R).
// Note that string constant template patterns are surrounded by "quotes"!
initialize() {
    defineConstant("s", {{s}});
    defineConstant("onset_time", {{onset_time}});
    defineConstant("target_pop", "{{target_pop}}");
    defineConstant("origin_pop", "{{origin_pop}}");

    // compose the path to a trajectory file based on given parameters
    defineConstant("traj_file",
                   PATH + "/" + "traj_" + target_pop + "_" + origin_pop + ".tsv");
}

function (void) add_mutation(void) {
    // sample one target carrier of the new mutation...
    target = sample(population(origin_pop).haplosomes, 1);
    // ... and add the mutation in the middle of it
    mut = target.addNewMutation(m1, s, position = asInteger(SEQUENCE_LENGTH / 2));

    // save the mutation for later reference
    defineGlobal("MUTATION", mut);

    write_log("adding beneficial mutation to population " + origin_pop);

    writeFile(traj_file, "tick\ttime\tfreq_origin\tfreq_target");
}

tick(onset_time) late() {
    // save simulation state in case we need to restart if the mutation is lost
    save_state();

    add_mutation();
}

tick(onset_time):SIMULATION_END late() {
    // the mutation is not segregating and is not fixed either -- we must restart
    if (!MUTATION.isSegregating & !MUTATION.isFixed) {
        write_log("mutation lost -- restarting");

        reset_state();

        add_mutation();
    }

    // compute the frequency of the mutation of interest and save it (if the
    // mutation is missing at this time, save its frequency as NA)
    freq_origin = "NA";
    freq_target = "NA";
    if (population(origin_pop, check = T))
      freq_origin = population(origin_pop).haplosomes.mutationFrequenciesInHaplosomes();
    if (population(target_pop, check = T))
      freq_target = population(target_pop).haplosomes.mutationFrequenciesInHaplosomes();

    writeFile(traj_file,
              community.tick + "\t" +
              model_time(community.tick) + "\t" +
              freq_origin + "\t" +
              freq_target, append = T);
}
)"

## -----------------------------------------------------------------------------
run_model <- function(origin_pop, onset_time) {
  extension <- substitute_values(
    extension_template, # <--- template SLiM code string directly (not as a file!)
    s = 0.1, onset_time = onset_time,
    origin_pop = origin_pop, target_pop = "EUR"
  )

  model <- compile_model(
    populations = list(afr, ooa, ehg, eur, ana, yam),
    gene_flow = gf, generation_time = 30,
    extension = extension
  )

  slim(model, sequence_length = 1e6, recombination_rate = 0,
       path = tempdir(), random_seed = 42)
}

run_model("EUR", onset_time = 15000)

head(load_traj("EUR"))

## -----------------------------------------------------------------------------
# create the ancestor of everyone and a chimpanzee outgroup
# (we set both N = 1 to reduce the computational time for this model)
chimp <- population("CH", time = 6.5e6, N = 1000)

# two populations of anatomically modern humans: Africans and Europeans
afr <- population("AFR", parent = chimp, time = 6e6, N = 10000)
eur <- population("EUR", parent = afr, time = 70e3, N = 5000)

# Neanderthal population splitting at 600 ky ago from modern humans
# (becomes extinct by 40 ky ago)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

# 5% Neanderthal introgression into Europeans between 55-50 ky ago
gf <- gene_flow(from = nea, to = eur, rate = 0.05, start = 55000, end = 45000)

model <- compile_model(
  populations = list(chimp, nea, afr, eur), gene_flow = gf,
  generation_time = 30
)

## ----introgression_model, echo=FALSE, fig.width=6, fig.height=4---------------
cowplot::plot_grid(
  plot_model(model, sizes = FALSE, proportions = TRUE),
  plot_model(model, sizes = FALSE, log = TRUE, proportions = TRUE),
  nrow = 1
)

## ----echo=FALSE---------------------------------------------------------------
eur <- population("EUR", time = 100e3, N = 10000)
nea <- population("NEA", time = 100e3, N = 1000, remove = 40000)

gf <- gene_flow(from = nea, to = eur, rate = 0.05, start = 55000, end = 54500)

## -----------------------------------------------------------------------------
extension <- r"(
initialize() {
  // model parameters to be substitute_values()'d from R below
  defineConstant("gene_length", {{gene_length}});
  defineConstant("n_genes", {{n_genes}});
  defineConstant("n_markers", {{n_markers}});
  defineConstant("introgression_time", {{introgression_time}});
  defineConstant("freq_file", PATH + "/{{freq_file}}");

  // total length of the genome to be simulated
  defineConstant("total_length", n_genes * gene_length);

  // positions of neutral Neanderthal markers along the genome
  defineConstant("neutral_pos", seq(0, total_length - 1, by = gene_length / n_markers));
  // positions of deleterious mutations in the center of each gene
  defineConstant("selected_pos", seq(gene_length / 2, total_length - 1, by = gene_length));
}

// Because we want to simulate non-neutral evolution, we have to provide a
// custom initialization callback -- slendr will use it to replace its default
// neutral genomic architecture (i.e. the initialize() {...} callback it uses
// by default for neutral simulations).
initialize() {
  initializeMutationType("m0", 0.5, "f", 0.0); // neutral Neanderthal marker mutation type
  initializeMutationType("m1", 0.5, "f", {{s}}); // deleterious Neanderthal mutation type
  initializeGenomicElementType("g1", m0, 1.0); // genomic type of 'genes'

  genes = c(); // gene start-end coordinates
  breakpoints = c(); // recombination breakpoints
  rates = c(); // recombination rates

  // compute coordinates of genes, as well as the recombination map breakpoints
  // between each gene
  start = 0;
  for (i in seqLen(n_genes)) {
    // end of the next gene
    end = start + gene_length - 1;
    genes = c(genes, start, end);

    // uniform recombination within a gene, followed by a 0.5 recombination breakpoint
    rates = c(rates, RECOMBINATION_RATE, 0.5);
    breakpoints = c(breakpoints, end, end + 1);

    // start of the following genes
    start = end + 1;
  }

  // odd elements --> starts of genes
  gene_starts = integerMod(seqAlong(genes), 2) == 0;
  // even elements --> ends of genes
  gene_ends = integerMod(seqAlong(genes), 2) == 1;

  // set up all the genes at once
  initializeGenomicElement(g1, genes[gene_starts], genes[gene_ends]);

  // set up the recombination map
  initializeRecombinationRate(rates, breakpoints);

  // no mutation rate (we will add neutral variation after the SLiM run)
  initializeMutationRate(0);
}

// Add Neanderthal-specific mutations one tick prior to the introgression
tick(introgression_time) - 1 late() {
  // get all Neanderthal chromosomes just prior to the introgression
  target = population("NEA").haplosomes;

  write_log("adding neutral Neanderthal markers");

  mutations = target.addNewDrawnMutation(m0, position = asInteger(neutral_pos));
  defineConstant("MARKERS", target.mutations);

  write_log("adding deleterious Neanderthal mutation");

  target.addNewDrawnMutation(m1, position = asInteger(selected_pos));
}

// At the end, write Neanderthal ancestry proportions along the genome to a file
// (in our R analysis code we will compare the results of this with the
// equivalent computation using a tree sequence)
SIMULATION_END late() {
  df = DataFrame(
    "gene", asInteger(MARKERS.position / gene_length),
    "pos", MARKERS.position,
    "freq", sim.mutationFrequencies(population("EUR"), MARKERS)
  );
  writeFile(freq_file, df.serialize(format = "tsv"));
}
)" %>%
  substitute_values(
    introgression_time = 55000, s = -0.003,
    gene_length = 5e6, n_genes = 200, n_markers = 100,
    freq_file = "frequencies.tsv"
  )

## -----------------------------------------------------------------------------
model <- compile_model(
  populations = list(nea, eur), gene_flow = gf,
  generation_time = 30,
  extension = extension
)

## -----------------------------------------------------------------------------
nea_samples <- schedule_sampling(model, times = 50000, list(nea, 1))
modern_samples <- schedule_sampling(model, times = 0, list(eur, 100))

samples <- rbind(nea_samples, modern_samples)

## ----echo=FALSE, eval = RERUN & run_vignette----------------------------------
# x = Sys.time()
# data_dir <- "~/Code/__archive/introgression_data/"
# slim(model, recombination_rate = 1e-8, samples = samples, path = data_dir)
# y = Sys.time()
# y - x # Time difference of 45.62365 mins

## ----eval=FALSE---------------------------------------------------------------
# slim(model, recombination_rate = 1e-8, samples = samples, path = "~/Code/__archive/introgression_data/")

## -----------------------------------------------------------------------------
n_genes <- 200
gene_length <- 5e6
window_length <- 100e3

## ----eval = run_vignette & file.exists("~/Code/__archive/introgression_data/frequencies.tsv")----
freqs <- read_tsv("~/Code/__archive/introgression_data/frequencies.tsv")

## ----introgression_freqs, eval = run_vignette & file.exists("~/Code/__archive/introgression_data/frequencies.tsv")----
freqs %>%
mutate(pos = pos %% 5e6) %>%
group_by(pos) %>%
summarise(freq = 100 * mean(freq)) %>%
ggplot() +
  geom_line(aes(pos, freq)) +
  geom_vline(xintercept = gene_length / 2, linetype = 2) +
  labs(x = "coordinate along a gene", y = "Neanderthal ancestry proportion [%]",
       title = "Proportion of Neanderthal ancestry in Europeans along 5Mb independent genes",
       subtitle = "(dashed line indicates introgressed deleterious Neanderthal allele)") +
  coord_cartesian(ylim = c(0, 3))

## ----eval = run_vignette & file.exists("~/Code/__archive/introgression_data/slim.trees")----
# load a tree sequence and extract the names of recorded individuals
ts <- ts_read(file = "~/Code/__archive/introgression_data/slim.trees", model)
samples <- ts_names(ts, split = "pop")

samples

# compute coordinates of sliding windows along the genome
windows <- seq(from = 0, to = n_genes * gene_length - 1, by = window_length)

head(windows)
tail(windows)

# compute divergence from the tree sequence in each window separately
divergence <- ts_divergence(ts, samples, windows = windows, mode = "branch")$divergence[[1]]

## ----eval = run_vignette & file.exists("~/Code/__archive/introgression_data/slim.trees")----
# compute average divergence at each position in a gene, and a 95% C.I.
div_df <- tibble(
  pop = "eur",
  pos = windows %% gene_length,
  div = divergence
) %>%
  group_by(pop, pos) %>%
  summarise(
    mean = mean(div), n = n(), std = sd(div),
    ci_low = mean - 2 * std / sqrt(n),
    ci_up = mean + 2 * std / sqrt(n)
  )

## ----introgression_divergence, eval = run_vignette & file.exists("~/Code/__archive/introgression_data/slim.trees")----
ggplot(div_df) +
  geom_ribbon(aes(pos, ymin = ci_low, ymax = ci_up), fill = "grey70") +
  geom_line(aes(pos, mean)) +
  geom_vline(aes(xintercept = gene_length / 2), linetype = 2) +
  labs(x = "coordinate along a gene", y = "divergence to Neanderthal",
       title = "Divergence of Europeans to a Neanderthal genome along 5Mb independent genes",
       subtitle = "(dashed line indicates introgressed deleterious Neanderthal allele)")

