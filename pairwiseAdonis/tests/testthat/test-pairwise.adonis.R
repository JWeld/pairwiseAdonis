data(iris)

test_that("significance codes follow the standard thresholds", {
  expect_equal(pairwiseAdonis:::sig_codes(c(0.2, 0.08, 0.03, 0.005, 0.0005, NA)),
               c("", ".", "*", "**", "***", ""))
})

test_that("basic call returns one row per pair with expected columns", {
  set.seed(1)
  res <- pairwise.adonis(iris[, 1:4], iris$Species, perm = 49)
  expect_s3_class(res, "pwadonis")
  expect_equal(nrow(res), 3)
  expect_named(res, c("pairs", "Df", "SumsOfSqs", "F.Model", "R2",
                      "p.value", "p.adjusted", "sig"))
  expect_equal(res$p.adjusted, p.adjust(res$p.value, "bonferroni"))
})

test_that("results agree with adonis2 on the same subset", {
  idx <- iris$Species %in% c("setosa", "versicolor")
  set.seed(7)
  res <- pairwise.adonis(iris[, 1:4], iris$Species, perm = 49)[1, ]
  set.seed(7)
  ad <- adonis2(vegdist(iris[idx, 1:4]) ~ Species, data = iris[idx, ], permutations = 49)
  expect_equal(res$SumsOfSqs, ad$SumOfSqs[1])
  expect_equal(res$F.Model, ad$F[1])
  expect_equal(res$R2, ad$R2[1])
  expect_equal(res$p.value, ad$`Pr(>F)`[1])
})

test_that("dist input and symmetric matrix input give identical results", {
  d <- vegdist(iris[, 1:4], "euclidean")
  set.seed(3); a <- pairwise.adonis(d, iris$Species, perm = 49)
  set.seed(3); b <- pairwise.adonis(as.matrix(d), iris$Species, perm = 49)
  set.seed(3); c <- pairwise.adonis(iris[, 1:4], iris$Species, sim.method = "euclidean", perm = 49)
  expect_equal(a, b)
  expect_equal(a, c)
})

test_that("reduce matches levels exactly, not as a regular expression", {
  set.seed(2)
  X <- matrix(abs(rnorm(90)), 30)
  f <- rep(c("a", "ab", "b"), each = 10)
  res <- pairwise.adonis(X, f, reduce = "a", perm = 19)
  expect_equal(res$pairs, c("a vs ab", "a vs b"))
  res2 <- pairwise.adonis(X, f, reduce = "a|b", perm = 19)
  expect_equal(nrow(res2), 3)
  expect_warning(pairwise.adonis(X, f, reduce = "a|zzz", perm = 19), "not found")
  expect_error(suppressWarnings(pairwise.adonis(X, f, reduce = "zzz", perm = 19)), "no pairwise comparison")
  # spaces around the pipes are ignored
  expect_equal(pairwise.adonis(X, f, reduce = " ab ", perm = 19)$pairs, c("a vs ab", "ab vs b"))
})

test_that("pairs follow the levels of a factor, otherwise the order of appearance", {
  set.seed(2)
  X <- matrix(abs(rnorm(90)), 30)
  f <- factor(rep(c("b", "a", "c"), each = 10), levels = c("a", "b", "c"))
  res <- pairwise.adonis(X, f, perm = 9)
  expect_equal(res$pairs, c("a vs b", "a vs c", "b vs c"))
  expect_equal(pairwise.adonis(X, as.character(f), perm = 9)$pairs, c("b vs a", "b vs c", "a vs c"))
  expect_type(res$pairs, "character")
  expect_type(res$sig, "character")
})

test_that("other p-value adjustments are applied", {
  res <- pairwise.adonis(iris[, 1:4], iris$Species, p.adjust.m = "holm", perm = 49)
  expect_equal(res$p.adjusted, p.adjust(res$p.value, "holm"))
  expect_error(pairwise.adonis(iris[, 1:4], iris$Species, p.adjust.m = "nope"), "should be one of")
})

test_that("daisy distances agree with daisy on the same subset", {
  idx <- iris$Species %in% c("setosa", "versicolor")
  set.seed(7)
  res <- pairwise.adonis(iris[, 1:4], iris$Species, sim.function = "daisy",
                         sim.method = "manhattan", perm = 49)[1, ]
  set.seed(7)
  ad <- adonis2(cluster::daisy(iris[idx, 1:4], metric = "manhattan") ~ Species,
                data = iris[idx, ], permutations = 49)
  expect_equal(res$F.Model, ad$F[1])
  expect_equal(res$p.value, ad$`Pr(>F)`[1])
  expect_error(pairwise.adonis(iris[, 1:4], iris$Species, sim.function = "daisy"), "must be one of")
  expect_error(pairwise.adonis(iris[, 1:4], iris$Species, sim.function = "nope"), "should be one of")
})

test_that("further arguments go to the distance function or to adonis2", {
  data(dune, dune.env, package = "vegan")
  idx <- dune.env$Management %in% c("BF", "HF")
  # binary = TRUE to vegdist: presence/absence Jaccard
  set.seed(9)
  res <- pairwise.adonis(dune, dune.env$Management, sim.method = "jaccard", binary = TRUE, perm = 49)[1, ]
  set.seed(9)
  ad <- adonis2(vegdist(dune[idx, ], "jaccard", binary = TRUE) ~ Management,
                data = dune.env[idx, ], permutations = 49)
  expect_equal(res$pairs, "BF vs HF")
  expect_equal(res$F.Model, ad$F[1])
  expect_equal(res$p.value, ad$`Pr(>F)`[1])
  # sqrt.dist = TRUE to adonis2
  set.seed(9)
  res <- pairwise.adonis(dune, dune.env$Management, sqrt.dist = TRUE, perm = 49)[1, ]
  set.seed(9)
  ad <- adonis2(vegdist(dune[idx, ]) ~ Management, data = dune.env[idx, ],
                permutations = 49, sqrt.dist = TRUE)
  expect_equal(res$F.Model, ad$F[1])
  expect_equal(res$p.value, ad$`Pr(>F)`[1])
  expect_error(pairwise.adonis(dune, dune.env$Management, permutations = 9), "perm")
  expect_error(pairwise.adonis(dune, dune.env$Management, "vegdist", "bray", "holm", NULL, 9, TRUE),
               "must be named")
})

test_that("input checks give informative errors", {
  expect_error(pairwise.adonis(iris[, 1:4], iris$Species[1:100]), "one entry per row")
  expect_error(pairwise.adonis(iris[1:50, 1:4], iris$Species[1:50]), "at least two levels")
  f <- iris$Species; f[1] <- NA
  expect_error(pairwise.adonis(iris[, 1:4], f), "missing values")
})

test_that("summary method is dispatched", {
  res <- pairwise.adonis(iris[, 1:4], iris$Species, perm = 9)
  expect_output(summary(res), "Result of pairwise.adonis")
})

test_that("a how() design with blocks is reduced to each pair", {
  set.seed(4)
  dat <- expand.grid(rep = gl(2, 1), NO3 = factor(c(0, 10, 30)), field = gl(3, 1))
  Y <- data.frame(A = with(dat, as.numeric(field) + as.numeric(NO3)) + abs(rnorm(18)),
                  B = with(dat, as.numeric(field) - as.numeric(NO3) + 6) + abs(rnorm(18)))
  h <- how(nperm = 49)
  setBlocks(h) <- dat$field
  set.seed(6); res <- pairwise.adonis(Y, dat$NO3, perm = h)
  idx <- dat$NO3 %in% c("0", "10")
  h2 <- how(nperm = 49)
  setBlocks(h2) <- dat$field[idx]
  set.seed(6); ad <- adonis2(vegdist(Y[idx, ]) ~ NO3, data = dat[idx, ], permutations = h2)
  expect_equal(res$F.Model[1], ad$F[1])
  expect_equal(res$p.value[1], ad$`Pr(>F)`[1])
  # designs that cannot be reduced are refused
  hp <- how(nperm = 49, plots = Plots(strata = dat$field))
  expect_error(pairwise.adonis(Y, dat$NO3, perm = hp), "plot-level strata")
  expect_error(pairwise.adonis(Y, dat$NO3, perm = shuffleSet(18, 9)), "permutation matrix")
  setBlocks(h) <- dat$field[1:6]
  expect_error(pairwise.adonis(Y, dat$NO3, perm = h), "one entry per observation")
  # restricted permutations within blocks and missing blocks are refused
  expect_error(pairwise.adonis(Y, dat$NO3, perm = how(within = Within(type = "series"), nperm = 9)),
               "restricted permutations")
  expect_error(pairwise.adonis(Y, dat$NO3,
                               perm = how(within = Within(type = "grid", nrow = 3, ncol = 6), nperm = 9)),
               "restricted permutations")
  b <- dat$field
  b[3] <- NA
  setBlocks(h) <- b
  expect_error(pairwise.adonis(Y, dat$NO3, perm = h), "missing values")
})

test_that("a block absent from a pair is handled", {
  set.seed(4)
  dat <- expand.grid(rep = gl(2, 1), NO3 = factor(c(0, 10, 30)), field = gl(3, 1))
  Y <- data.frame(A = with(dat, as.numeric(field) + as.numeric(NO3)) + abs(rnorm(18)),
                  B = with(dat, as.numeric(field) - as.numeric(NO3) + 6) + abs(rnorm(18)))
  # block "4" contains only observations of level 30
  blk <- factor(ifelse(dat$NO3 == "30" & dat$field == "3", "4", as.character(dat$field)))
  h <- how(nperm = 49)
  setBlocks(h) <- blk
  set.seed(6); res <- pairwise.adonis(Y, dat$NO3, perm = h)
  idx <- dat$NO3 %in% c("0", "10")
  h2 <- how(nperm = 49)
  setBlocks(h2) <- blk[idx]
  set.seed(6); ad <- adonis2(vegdist(Y[idx, ]) ~ NO3, data = dat[idx, ], permutations = h2)
  expect_equal(res$F.Model[1], ad$F[1])
  expect_equal(res$p.value[1], ad$`Pr(>F)`[1])
})
