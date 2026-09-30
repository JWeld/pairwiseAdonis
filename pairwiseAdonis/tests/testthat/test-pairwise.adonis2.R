set.seed(11)
dat <- expand.grid(rep = gl(2, 1), NO3 = factor(c(0, 10, 30)), field = gl(3, 1))
dat$NO3num <- as.numeric(as.character(dat$NO3))
Y <- data.frame(A = with(dat, as.numeric(field) + as.numeric(NO3) + 2) + abs(rnorm(18)) / 2,
                B = with(dat, as.numeric(field) - as.numeric(NO3) + 6) + abs(rnorm(18)) / 2)

test_that("basic call returns parent call and one anova table per pair", {
  res <- pairwise.adonis2(Y ~ NO3, data = dat, nperm = 49)
  expect_s3_class(res, "pwadstrata")
  expect_named(res, c("parent_call", "0_vs_10", "0_vs_30", "10_vs_30"))
  expect_s3_class(res[["0_vs_10"]], "anova")
  expect_true(any(grepl("Number of permutations: 49", attr(res[["0_vs_10"]], "heading"))))
})

test_that("nperm is respected with and without strata", {
  res <- pairwise.adonis2(Y ~ NO3, data = dat, nperm = 29)
  expect_true(any(grepl("Number of permutations: 29", attr(res[[2]], "heading"))))
  res <- pairwise.adonis2(Y ~ NO3, data = dat, strata = "field", nperm = 29)
  expect_true(any(grepl("Number of permutations: 29", attr(res[[2]], "heading"))))
  expect_true(any(grepl("Blocks", attr(res[[2]], "heading"))))
})

test_that("results with strata equal a manual adonis2 call on the subset", {
  idx <- dat$NO3 %in% c("0", "10")
  perm <- how(nperm = 99)
  setBlocks(perm) <- dat$field[idx]
  set.seed(5); a <- pairwise.adonis2(Y ~ NO3, data = dat, strata = "field", nperm = 99)[["0_vs_10"]]
  set.seed(5); b <- adonis2(Y[idx, ] ~ NO3, data = dat[idx, ], permutations = perm)
  expect_equal(as.data.frame(a), as.data.frame(b), ignore_attr = TRUE)
})

test_that("dist and symmetric matrix responses give identical results", {
  d <- vegdist(Y, "euclidean")
  set.seed(3); a <- pairwise.adonis2(d ~ NO3, data = dat, nperm = 49)
  set.seed(3); b <- pairwise.adonis2(as.matrix(d) ~ NO3, data = dat, nperm = 49)
  set.seed(3); c <- pairwise.adonis2(Y ~ NO3, data = dat, nperm = 49, method = "euclidean")
  expect_equal(as.data.frame(a[[2]]), as.data.frame(b[[2]]), ignore_attr = TRUE)
  expect_equal(as.data.frame(a[[2]]), as.data.frame(c[[2]]), ignore_attr = TRUE)
})

test_that("function calls on the right hand side are preserved", {
  res <- pairwise.adonis2(Y ~ factor(NO3num), data = dat, nperm = 19, by = "terms")
  expect_named(res, c("parent_call", "0_vs_10", "0_vs_30", "10_vs_30"))
  expect_equal(rownames(res[[2]])[1], "factor(NO3num)")
})

test_that("extra terms are kept and by = 'terms' is forwarded", {
  res <- pairwise.adonis2(Y ~ NO3/field, data = dat, strata = "field", nperm = 19, by = "terms")
  expect_equal(rownames(res[[2]]), c("NO3", "NO3:field", "Residual", "Total"))
})

test_that("input checks give informative errors", {
  expect_error(pairwise.adonis2(Y ~ NO3, data = dat[1:17, ]), "same number of rows")
  expect_error(pairwise.adonis2(Y ~ NO3, data = dat, strata = "nope"), "column of 'data'")
  expect_error(pairwise.adonis2(Y ~ NO3, data = dat, permutations = 9), "nperm")
  datNA <- dat; datNA$NO3[1] <- NA
  expect_error(pairwise.adonis2(Y ~ NO3, data = datNA), "missing values in the grouping")
  expect_error(pairwise.adonis2(Y[dat$NO3 == "0", ] ~ NO3, data = dat[dat$NO3 == "0", ]), "at least two levels")
})

test_that("summary method is dispatched", {
  res <- pairwise.adonis2(Y ~ NO3, data = dat, nperm = 9)
  expect_output(summary(res), "Result of pairwise.adonis2")
})

test_that("functions defined by the user can be used on the right hand side", {
  # a function visible from the global environment, as in an interactive session
  assign("my_tr", function(z) factor(z), envir = globalenv())
  on.exit(rm("my_tr", envir = globalenv()))
  res <- pairwise.adonis2(Y ~ my_tr(NO3num), data = dat, nperm = 9, by = "terms")
  expect_named(res, c("parent_call", "0_vs_10", "0_vs_30", "10_vs_30"))
  expect_equal(rownames(res[[2]])[1], "my_tr(NO3num)")
})

test_that("objects local to the caller are found when adonis2 itself supports it", {
  # adonis2() replaces the formula environment up to vegan 2.7 and keeps it
  # from vegan 2.8-0; the wrapper must not be more restrictive than adonis2()
  adonis2_keeps_env <- function() {
    h <- function() {
      loc_tr <- function(z) factor(z)
      adonis2(Y ~ loc_tr(NO3num), data = dat, permutations = 9)
    }
    !inherits(try(h(), silent = TRUE), "try-error")
  }
  skip_if_not(adonis2_keeps_env(), "adonis2() in this vegan version drops the formula environment")
  f <- function() {
    loc_tr <- function(z) factor(z)
    pairwise.adonis2(Y ~ loc_tr(NO3num), data = dat, nperm = 9, by = "terms")
  }
  expect_equal(rownames(f()[[2]])[1], "loc_tr(NO3num)")
  g <- function() {
    cutoff <- 1.5
    pairwise.adonis2(Y ~ NO3 + I(as.numeric(rep) > cutoff), data = dat, nperm = 9, by = "terms")
  }
  expect_equal(rownames(g()[["0_vs_30"]])[1:2], c("NO3", "I(as.numeric(rep) > cutoff)"))
})

test_that("NA in a secondary right hand side term keeps response and data aligned", {
  datNA <- dat
  datNA$field[18] <- NA   # an observation of level 30
  # pair not involving the NA row: identical to the complete data
  set.seed(8); a <- pairwise.adonis2(Y ~ NO3 + field, data = datNA, nperm = 49,
                                     na.action = na.omit)[["0_vs_10"]]
  set.seed(8); b <- pairwise.adonis2(Y ~ NO3 + field, data = dat, nperm = 49)[["0_vs_10"]]
  expect_equal(as.data.frame(a), as.data.frame(b), ignore_attr = TRUE)
  # pair involving the NA row: identical to adonis2 with na.omit on the manual subset
  idx <- dat$NO3 %in% c("0", "30")
  # (this pair is not the first computed, so the permutation p-values use a
  # different random stream; compare the deterministic columns)
  a <- pairwise.adonis2(Y ~ NO3 + field, data = datNA, nperm = 49, na.action = na.omit)[["0_vs_30"]]
  b <- adonis2(Y[idx, ] ~ NO3 + field, data = datNA[idx, ], permutations = 49, na.action = na.omit)
  expect_equal(as.data.frame(a)[, 1:4], as.data.frame(b)[, 1:4], ignore_attr = TRUE)
  expect_equal(a$Df[length(a$Df)], sum(idx) - 2)
  # same with a distance matrix response
  d <- vegdist(Y, "euclidean")
  a <- pairwise.adonis2(d ~ NO3 + field, data = datNA, nperm = 49, na.action = na.omit)[["0_vs_30"]]
  b <- adonis2(as.dist(as.matrix(d)[idx, idx]) ~ NO3 + field, data = datNA[idx, ],
               permutations = 49, na.action = na.omit)
  expect_equal(as.data.frame(a)[, 1:4], as.data.frame(b)[, 1:4], ignore_attr = TRUE)
  # default na.fail still errors for the affected pair only
  expect_error(pairwise.adonis2(Y ~ NO3 + field, data = datNA, nperm = 9), "missing values")
})
