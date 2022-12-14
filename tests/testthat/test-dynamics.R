map <- readRDS("map.rds")

test_that("move requires overlap between 0 and 1", {
  p <- population(name = "pop2", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, overlap = -1),
               "The required overlap between subsequent spatial maps must be a number between 0 and 1", fixed = TRUE)
  expect_error(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, overlap = 0),
               "The required overlap between subsequent spatial maps must be a number between 0 and 1", fixed = TRUE)
  expect_silent(suppressMessages(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, overlap = 0.3)))
})

test_that("move requires a positive snapshot number", {
  p <- population(name = "pop2", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, snapshots = -1),
               "The number of snapshots must be a non-negative integer", fixed = TRUE)
  expect_error(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, snapshots = 0),
               "The number of snapshots must be a non-negative integer", fixed = TRUE)
  expect_silent(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, snapshots = 3, verbose = FALSE))
})

test_that("move generates the correct number of snapshots", {
  p <- population(name = "pop2", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  nsnapshots <- 3
  moved_p <- move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, snapshots = nsnapshots, verbose = FALSE)
  expect_true(nrow(moved_p) == nsnapshots + 2)
})

test_that("only positive population sizes allowed", {
  p <- population(name = "pop", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(resize(p, N = -1, time = 1000, how = "step"),
               "resize\\(\\): Only positive, non-zero population sizes are allowed")
})

test_that("only three values of population size changes allowed", {
  p <- population(name = "pop", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(resize(p, N = 10, how = "asdf", time = 1000),
               "resize\\(\\): Only 'step' or 'exponential' are allowed")
})
