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
