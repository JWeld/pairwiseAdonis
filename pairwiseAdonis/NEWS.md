# pairwiseAdonis 0.5.0

Further fixes after a review of 0.4.2. Tested against vegan 2.6-4, and against
the `adonis2()` code of vegan 2.7-6 and of the development version 2.8-0.

This release follows 0.4.2 of this fork. The version jumps to 0.5.0 because
default results of `pairwise.adonis2()` change (adjusted p-values), and to
leave the 0.4.x numbers to the original package
(https://github.com/pmartinezarbizu/pairwiseAdonis).

## Changes in results

* `pairwise.adonis2()` now adjusts p-values for multiple comparisons, as
  `pairwise.adonis()` does. New argument `p.adjust.m` (default
  `'bonferroni'`): the p-value in the first row of the table of each pair is
  adjusted for the number of pairs and added as column `Pr(adj)`; the
  significance codes then refer to the adjusted p-values and the unadjusted
  ones stay in `Pr(>F)`. Use `p.adjust.m = 'none'` for the tables of
  earlier versions. The first row is the test of the grouping variable when
  it is the only term or with `by = "terms"`; with several terms and
  `by = NULL` it is the omnibus test of the whole model.
* Both functions: when the grouping variable is a factor, the pairs follow
  the order of its levels instead of the order in which the levels first
  appear in the data. Pairs are therefore computed in a different order for
  such data, and results obtained with `set.seed()` change accordingly.

## Bug fixes

* `pairwise.adonis2()`: `na.action = na.omit` together with `strata` failed
  ("Number of observations and length of Block 'strata' do not match").
  Observations with missing values on the right-hand side are now dropped
  before the pairs are formed when `na.action` allows it, so that the
  response, the data and the strata stay aligned. Missing values in the
  grouping variable are dropped too with `na.action = na.omit`; with the
  default `na.fail` the error message names the problem.
* Both functions: missing values in `strata` (or in the blocks of a `how()`
  design) gave a cryptic "subscript out of bounds" error. They now give an
  informative error.
* `pairwise.adonis2()`: a numeric first variable on the right-hand side was
  silently treated as a grouping variable, giving one comparison for every
  pair of distinct values (e.g. 78 comparisons for a pH variable). This now
  gives a warning that names the variable and the number of comparisons;
  use `factor(x)` in the formula for numeric group codes. A matrix-valued
  first variable (e.g. `poly(x, 2)`) gives an error.
* `pairwise.adonis()`: `how()` designs with restricted permutations within
  blocks (`Within(type = "series")` or `"grid"`) were applied to the subset
  of each pair, where their spatial or temporal structure no longer holds
  (series ran silently, grid failed with an unhelpful error). They are now
  refused with an informative error, like plot-level designs.

## vegan 2.8-0

* `pairwise.adonis2()` accepts `Condition()` terms (partial models, new in
  vegan 2.8-0's `adonis2()`). They are passed to `adonis2()` and never used
  as grouping variable.

## Other

* `pairwise.adonis()` gained `...`: `sqrt.dist`, `add` and `parallel` are
  passed to `adonis2()`, all other arguments to the distance function, for
  example `binary = TRUE` for presence/absence indices with `vegdist()`.
* `pairwise.adonis()`: spaces around the pipes in `reduce` are ignored;
  `sim.function` and `p.adjust.m` are checked before any comparison is
  computed; `sim.function = 'daisy'` with a method `daisy()` does not
  support (such as the default `'bray'`) gives an informative error; the
  `pairs` and `sig` columns are character also in R < 4.0.
* `pairwise.adonis2()`: informative errors when `x` is not a formula with a
  response, when `data` is not a data frame, or when there is no grouping
  variable. The `adonis2()` tables show the pair and the strata column in
  their heading instead of internal object names, the parent call records
  further arguments such as `by = "terms"`, and the result has a `print()`
  method.
* `DESCRIPTION`: Sylvain Monteux is listed as author, and the maintainer of
  this fork as maintainer.

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
* `pairwise.adonis2()`: the environment of the model formula is kept, so
  functions and objects of the caller used on the right-hand side are found
  wherever `adonis2()` itself finds them (from vegan 2.8-0 `adonis2()` keeps
  the formula environment too).
* `pairwise.adonis()`: a permutation design from `permute::how()` with
  blocks is reduced to the observations of each pair. Designs with
  plot-level strata and permutation matrices are refused with an
  informative error, because each pair uses a different subset of the
  observations.
* `summary()` methods for both result classes are now registered and work.
* Informative errors when `factors`/`data` do not match the number of
  observations, when there are fewer than two levels, or when the grouping
  variable contains `NA`.

## Other

* `DESCRIPTION`: `cluster`, `stats` and `utils` moved to `Imports`; `vegan`
  and `permute` remain in `Depends` so that `vegdist()`, `how()` and
  `setBlocks()` are available after `library(pairwiseAdonis)`.
* Added a testthat test suite.
