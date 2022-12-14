skip_if(!slendr:::check_env_present())

set.seed(42)

script_file <- tempfile()
ts_file <- tempfile()
loc_file <- tempfile()

writeLines(sprintf('
initialize() {
  initializeSLiMOptions(keepPedigrees = T, dimensionality = "xy");
  initializeTreeSeq();
  initializeMutationRate(0);
  initializeMutationType("m1", 0.5, "f", 0.0);
  initializeGenomicElementType("g1", m1, 1.0);
  initializeGenomicElement(g1, 0, 1e6);
  initializeRecombinationRate(1e-8);
}
1 late() {
  sim.addSubpop("p1", 100);
  p1.individuals.x = runif(p1.individualCount);
  p1.individuals.y = runif(p1.individualCount);
}
modifyChild() {
  do child.x = parent1.x + rnorm(1, 0, 0.02);
  while ((child.x < 0.0) | (child.x > 1.0));

  do child.y = parent1.y + rnorm(1, 0, 0.02);
  while ((child.y < 0.0) | (child.y > 1.0));

  return T;
}
5000 late() {
  sim.treeSeqOutput("%s");
  for (ind in sim.subpopulations.individuals) {
    writeFile("%s", paste(ind.spatialPosition, ind.pedigreeID, sep = "\t"), append = T);
  }
}
', ts_file, loc_file), script_file)

system2("slim", script_file, stdout = FALSE)

suppressMessages(ts <- ts_load(ts_file, simplify = TRUE))

test_that("non-slendr SLiM tree sequence locations are correctly loaded", {
  data <- ts_nodes(ts, sf = FALSE) %>%
    dplyr::arrange(pedigree_id) %>%
    dplyr::select(x, y, pedigree_id) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(pedigree_id)) %>%
    as.data.frame()
  locations <- readr::read_tsv(
    loc_file, col_types = "ddi",
    col_names = c("x", "y", "pedigree_id")
  ) %>% dplyr::arrange(pedigree_id)

  expect_true(all.equal(data$x, locations$x, tolerance = 0.00001))
  expect_true(all.equal(data$y, locations$y, tolerance = 0.00001))
  expect_true(all(data$pedigree_id == locations$pedigree_id))
})

