## ---- include = FALSE---------------------------------------------------------
env_present <- slendr:::is_slendr_env_present()

knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  dpi = 80,
  eval = Sys.which("slim") != "" && env_present && Sys.getenv("RUNNER_OS") != "macOS"
)

## ---- collapse = TRUE, message = FALSE----------------------------------------
library(slendr)

library(dplyr)
library(ggplot2)

init_env()

seed <- 314159
set.seed(seed)

## ---- results = FALSE---------------------------------------------------------
# simulated world map
map <- world(
  xrange = c(-13, 70), # min-max longitude
  yrange = c(18, 65),  # min-max latitude
  crs = "EPSG:3035"    # coordinate reference system (CRS) for West Eurasia
)

# couple of broad geographic regions
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

# define population histories

# African ancestral population
afr <- population(
  "AFR", time = 52000, N = 3000,
  map = map, polygon = africa
)

# population of the first migrants out of Africa
ooa <- population(
  "OOA", parent = afr, time = 51000, N = 500, remove = 25000,
  center = c(33, 30), radius = 400e3
) %>%
  move(
    trajectory = list(c(40, 30), c(50, 30), c(60, 40)),
    start = 50000, end = 40000, snapshots = 20
  )

# Eastern hunter-gatherers
ehg <- population(
  "EHG", parent = ooa, time = 28000, N = 1000, remove = 6000,
  polygon = list(
    c(26, 55), c(38, 53), c(48, 53), c(60, 53),
    c(60, 60), c(48, 63), c(38, 63), c(26, 60))
)

# European population
eur <- population(name = "EUR", parent = ehg, time = 25000, N = 2000, polygon = europe)

# Anatolian farmers
ana <- population(
  name = "ANA", time = 28000, N = 3000, parent = ooa, remove = 4000,
  center = c(34, 38), radius = 500e3, polygon = anatolia
) %>%
  expand_range(
    by = 2500e3, start = 10000, end = 7000,
    polygon = join(europe, anatolia), snapshots = 20
  ) # expand the range by 2.500 km

# Yamnaya steppe population
yam <- population(
  name = "YAM", time = 7000, N = 500, parent = ehg, remove = 2500,
  polygon = list(c(26, 50), c(38, 49), c(48, 50),
                 c(48, 56), c(38, 59), c(26, 56))
) %>%
  move(trajectory = list(c(15, 50)), start = 5000, end = 3000, snapshots = 10)

# geneflow events
gf <- list(
  gene_flow(from = ana, to = yam, rate = 0.5, start = 6000, end = 5000, overlap = FALSE),
  gene_flow(from = ana, to = eur, rate = 0.5, start = 8000, end = 6000),
  gene_flow(from = yam, to = eur, rate = 0.75, start = 4000, end = 3000)
)

# compile the spatial model
model <- compile_model(
  populations = list(afr, ooa, ehg, eur, ana, yam),
  gene_flow = gf,
  generation_time = 30, resolution = 10e3,
  competition = 150e3, mating = 120e3, dispersal = 90e3
)

## ---- demographic_model-------------------------------------------------------
plot_model(model)

## ---- spatial_maps------------------------------------------------------------
plot_map(afr, ooa, ehg, eur, ana, yam)

## -----------------------------------------------------------------------------
# one ancient individual every two thousand years
ancient <- schedule_sampling(model,
                    times = seq(40000, 1, by = -500),
                    list(ooa, 1), list(ehg, 1), list(eur, 1),
                    list(ana, 1), list(yam, 1))

# present-day Africans and Europeans
present <- schedule_sampling(model, times = 0, list(afr, 5), list(eur, 30))

samples <- rbind(ancient, present)

## -----------------------------------------------------------------------------
ts <- slim(
  model, sequence_length = 100e3, recombination_rate = 1e-8, burnin = 200e3,
  samples = samples, method = "batch", random_seed = 314159, max_attempts = 1
) %>%
  ts_recapitate(recombination_rate = 1e-8, Ne = 10000, random_seed = seed) %>%
  ts_simplify()

ts

## -----------------------------------------------------------------------------
data <- ts_nodes(ts)

## -----------------------------------------------------------------------------
class(data)

## -----------------------------------------------------------------------------
data

## -----------------------------------------------------------------------------
map

## ---- slendr_map_ts-----------------------------------------------------------
sampled_data <- ts_nodes(ts) %>% filter(sampled)

ggplot() +
  geom_sf(data = map, fill = "lightgray", color = NA) +
  geom_sf(data = sampled_data, aes(shape = pop, color = time)) +
  ggtitle("Locations of simulated sampled individuals") +
  scale_color_continuous(type = "viridis") +
  theme_bw()

## -----------------------------------------------------------------------------
epochs <- sampled_data %>%
  mutate(epoch = cut(time, breaks = c(40000, 30000, 10000, 4000, 0)),
         epoch = ifelse(is.na(epoch), 0, epoch),
         epoch = factor(epoch, labels = c("present", "(present, 4 ky]", "(4 ky, 10 ky]",
                                          "(10 ky, 30 y]", "(30 ky, 40 ky]")))

## ---- map_epochs--------------------------------------------------------------
ggplot() +
  geom_sf(data = map, fill = "lightgray", color = NA) +
  geom_sf(data = epochs, aes(shape = pop, color = pop)) +
  facet_wrap(~ epoch) +
  ggtitle("Locations of simulated sampled individuals in different epochs") +
  theme_bw()

## -----------------------------------------------------------------------------
ind <- "EUR_67"

lineages <- ts_ancestors(ts, ind, verbose = TRUE)

## -----------------------------------------------------------------------------
lineages

## -----------------------------------------------------------------------------
filter(lineages, level == 1)

## ---- include = FALSE---------------------------------------------------------
counts <- filter(lineages, name == ind, level == 1) %>% as.data.frame() %>% count(node_id)

## ----level1_branches----------------------------------------------------------
level1_branches <- ts_ancestors(ts, "EUR_67") %>% filter(level == 1)

ggplot() +
  geom_sf(data = map, fill = "lightgray", color = NA) +
  geom_sf(data = level1_branches[, ]$child_location, shape = 13, size = 3, color = "red") +
  geom_sf(data = level1_branches[, ]$connection, linetype = 3) +
  geom_sf(data = level1_branches[, ]$parent_location, shape = 20, color = "blue") +
  theme_bw() +
  ggtitle("Parent nodes (blue) of a focal individual (red)")

## -----------------------------------------------------------------------------
as_tibble(level1_branches)[, c("name", "node_id", "child_id", "parent_id", "left_pos", "right_pos")]

## ---- eval = FALSE, echo = FALSE----------------------------------------------
#  # pdf("~/Desktop/x.pdf")
#  # for (i in ts_samples(ts) %>% filter(pop == "EUR") %>% pull(name)) {
#  
#  # ancestors <- ts_ancestors(ts, i, verbose = TRUE)
#  #   p_ancestors <- ggplot() +
#  #   geom_sf(data = map) +
#  #   geom_sf(data = ancestors, size = 0.5, aes(alpha = parent_time)) +
#  #   geom_sf(data = sf::st_set_geometry(ancestors, "parent_location"),
#  #           aes(shape = parent_pop, color = parent_pop)) +
#  #   geom_sf(data = filter(ts_nodes(ts), name == i), size = 3) +
#  #   coord_sf(expand = 0) +
#  #   labs(x = "longitude", y = "latitude") +
#  #   theme_bw() +
#  #   facet_grid(. ~ node_id) +
#  #   ggtitle(i) +
#  #   theme(legend.position = "none"); print(p_ancestors)
#  
#  # }
#  # dev.off()

## ----plot_ancestors_time------------------------------------------------------
ggplot() +
  geom_sf(data = map) +
  geom_sf(data = lineages, size = 0.5, alpha = 0.2) +
  geom_sf(data = sf::st_set_geometry(lineages, "parent_location"),
          aes(shape = parent_pop, color = parent_pop)) +
  geom_sf(data = filter(ts_nodes(ts), name == ind), size = 3) +
  guides(alpha = "none") +
  coord_sf(expand = 0) +
  labs(x = "longitude", y = "latitude") +
  facet_grid(. ~ node_id) +
  ggtitle("Ancestry encoded by two nodes (chromosomes) of EUR_67")

## ---- include = F-------------------------------------------------------------
nodes <- unique(lineages[!is.na(lineages$name), ]$node_id)

## -----------------------------------------------------------------------------
ts_samples(ts) %>% filter(pop == "ANA")

## ---- eval = FALSE, echo = FALSE----------------------------------------------
#  # pdf("~/Desktop/z.pdf")
#  # for (i in ts_samples(ts) %>% filter(pop == "ANA") %>% pull(name)) {
#  
#  # ancestors <- ts_ancestors(ts, i, verbose = TRUE)
#  #   p_ancestors <- ggplot() +
#  #   geom_sf(data = map) +
#  #   geom_sf(data = ancestors, size = 0.5, aes(alpha = parent_time)) +
#  #   geom_sf(data = sf::st_set_geometry(ancestors, "parent_location"),
#  #           aes(shape = parent_pop, color = parent_pop)) +
#  #   geom_sf(data = filter(ts_nodes(ts), name == i), size = 3) +
#  #   coord_sf(expand = 0) +
#  #   labs(x = "longitude", y = "latitude") +
#  #   theme_bw() +
#  #   facet_grid(. ~ node_id) +
#  #   ggtitle(i) +
#  #   theme(legend.position = "none"); print(p_ancestors)
#  
#  # }
#  # dev.off()

## ---- include = FALSE---------------------------------------------------------
ana_ind <- ts_samples(ts) %>% filter(name == "ANA_45")

## ---- plot_ancestors_levels_ana-----------------------------------------------
lineages <- ts_ancestors(ts, "ANA_45")

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = lineages, size = 0.5, alpha = 0.2) +
  geom_sf(data = sf::st_set_geometry(lineages, "parent_location"),
          aes(shape = parent_pop, color = parent_pop)) +
  geom_sf(data = filter(ts_nodes(ts), name == "ANA_45"), size = 3) +
  guides(alpha = "none") +
  coord_sf(expand = 0) +
  labs(x = "longitude", y = "latitude") +
  facet_grid(. ~ node_id) +
  ggtitle("Ancestry encoded by two nodes (chromosomes) of ANA_45")

## -----------------------------------------------------------------------------
lineages <-
  ts_samples(ts) %>%
  pull(name) %>%
  ts_ancestors(ts, x = .)

select(lineages, connection, child_time, parent_time)

## -----------------------------------------------------------------------------
distances <- lineages %>%
  mutate(branch_length = abs(parent_time - child_time) / model$generation_time,
         distance = sf::st_length(connection) %>% units::set_units(km) %>% as.numeric(),
         speed = distance / branch_length,
         epoch = cut(parent_time, breaks = c(Inf, seq(60000, 0, by = -3000)), dig.lab = 10, include.lowest = TRUE)) %>%
  as_tibble() %>% # strip away the spatial annotation
  select(name, pop, node_id, branch_length, distance, speed, parent_pop, parent_time, child_pop, child_time, epoch)

## -----------------------------------------------------------------------------
distances_long <- distances %>%
  filter(child_time < 60000) %>%
  filter(!pop %in% c("AFR", "OOA")) %>%
  tidyr::pivot_longer(cols = c(distance, speed),
                      names_to = "stat",
                      values_to = "value") %>%
  mutate(facet = case_when(
    stat == "distance" ~ "absolute distance of a node from parent",
    stat == "speed" ~ "distance traveled by a node per generation"))

## ---- smooth_distance_fits----------------------------------------------------
distances_long %>%
  ggplot(aes(child_time, value, color = child_pop)) +
  geom_smooth(method = "loess", aes(group = child_pop)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5) +
  labs(y = "kilometers", x = "time [years ago]") +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        legend.position = "bottom") +
  facet_wrap(~ facet, scales = "free_y") +
  guides(color = guide_legend("ancestral node population"))

