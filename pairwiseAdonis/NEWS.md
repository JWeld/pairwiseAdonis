# pairwiseAdonis 0.4.2

Compatibility with vegan 2.7 and bug fixes. Checked against vegan 2.6-4 and
vegan 2.7-6 (current CRAN release).

## vegan compatibility

* `adonis()` is defunct in vegan >= 2.7-1. The package no longer imports it;
  both functions have used `adonis2()` since version 0.4.
* `adonis2()` defaults to `by = NULL` (an omnibus test of the whole model)
  since vegan 2.6-8. With a single grouping factor the pairwise results are
  unchanged. For models with several terms (e.g. `Y ~ NO3/field`) pass
  `by = "terms"` to `pairwise.adonis2()` to get one row per term as before.
  This is now documented.
* `pairwise.adonis2()` results now keep the full `adonis2()` anova table,
  including its heading (call, number of permutations and blocks), instead
  of the first five columns only.

## Bug fixes

* `pairwise.adonis()`: the significance codes in the `sig` column were one
  step too strict (e.g. `'*'` was only given for adjusted p <= 0.01) unless
  `reduce` was used. They now follow the standard codes
  `0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1` in all cases.
* `pairwise.adonis()`: `reduce` matched level names as a regular expression,
  so `reduce = 'a'` also kept pairs of a level called `ab`. Levels are now
  matched exactly (still separated by `|`). Only the requested pairs are
  computed.
* `pairwise.adonis2()`: `nperm` was ignored when `strata` was not given
  (999 permutations were always used).
* Both functions: a symmetric square matrix (e.g. `as.matrix(vegdist(x))`)
  was treated as a community table and re-transformed with Bray-Curtis,
  giving wrong results silently. It is now treated as a distance matrix,
  as `adonis2()` does.
* `pairwise.adonis2()`: the right-hand side of the formula was rebuilt by
  pasting its parts, which lost function calls such as `~ factor(x)`.
  The original right-hand side is now used as is.
* `pairwise.adonis2()`: rows with missing values in the right-hand side
  variables were dropped when building the model frame, which could
  misalign the response and the data. Rows stay aligned; missing values
  are handled by `adonis2()` through its `na.action` argument.
* `summary()` methods for both result classes are now registered and work.
* Informative errors when `factors`/`data` do not match the number of
  observations, when there are fewer than two levels, or when the grouping
  variable contains `NA`.

## Other

* `DESCRIPTION`: `cluster`, `stats` and `utils` moved to `Imports`; `vegan`
  and `permute` remain in `Depends` so that `vegdist()`, `how()` and
  `setBlocks()` are available after `library(pairwiseAdonis)`.
* Added a testthat test suite.
