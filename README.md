# pairwiseAdonis
# version 0.4.2 includes 2 functions
pairwise.adonis

pairwise.adonis2

# pairwise.adonis
This is a wrapper function for multilevel pairwise comparison using adonis2 (~Permanova) from package 'vegan'. The function returns adjusted p-values using p.adjust(). It does not accept interaction between factors neither strata.

# pairwise.adonis2
This function accepts strata

NOTE: This is still a developing version -- Please validate your results.
I would appreciate feed back.

update 28 April 2020:
Function now adapts the permutation matrix for each combination of factors before applying adonis.
The p-value when using strata (block) looks correct now. Check syntax in example below.

This function accepts a model formula like in adonis from vegan. You can use interactions between factors and define strata to constrain permutations. For pairwise comparison a list of unique pairwise combination of factors is produced. Then for each pair, following objects are reduced accordingly to include only the subset of cases belonging to the pair:

- the left hand side of the formula (dissimilarity matrix or community matrix)

- the right hand side of the formula (factors)

- the strata if used.

The reduced data are passed to adonis and the summary of the anova table for each pair is saved in a list, together with the anova table of the full model and the original 'parent call'.

update 30.09.2026 (version 0.4.2)
Checked against vegan 2.6-4 and vegan 2.7-6 (adonis() is defunct since vegan 2.7-1; the package only uses adonis2()).
Bug fixes: significance codes in pairwise.adonis were too strict, `reduce` matched level names as a regular expression,
`nperm` was ignored in pairwise.adonis2 without strata, and symmetric distance matrices passed as plain matrices were
silently re-transformed with Bray-Curtis. See NEWS.md in the package for the full list.

NOTE on vegan >= 2.6-8: adonis2 now defaults to `by = NULL` (one omnibus test of the whole model per pair).
For models with several terms, e.g. `Y ~ NO3/field`, pass `by = "terms"` to pairwise.adonis2 to get one row per term:

```pairwise.adonis2(Y ~ NO3/field, data = dat, strata = 'field', by = "terms")```

update 30.08.21
Both functions now use adonis2 instead of adonis. This will solve some problems when loading DescTools. The functions should now work with Phyloseq objects.

Thanks to @lkoest12 for raising this and @JFMSilva for the solution 
https://github.com/joey711/phyloseq/issues/1457#issuecomment-880093708




Here a sample code provided by Andrea Pietruska for phyloseq:

For Jaccard and Bray Curtis:

```ps_tr <- microbiome::transform(phyloseq_object, "clr")```

```ps_dist_matrix <- phyloseq::distance(phyloseq_object, method ="jaccard")```

```vegan::adonis2(ps_dist_matrix ~ phyloseq::sample_data(ps_tr)$treatment)```

```pairwise.adonis(ps_dist_matrix, phyloseq::sample_data(ps_tr)$treatment)```


For Unifrac:

```weighted_unifrac = UniFrac(physeq=phyloseq_object, weighted=T, normalized=T, parallel=T, fast=T)```

```vegan::adonis2(weighted_unifrac ~ phyloseq::sample_data(ps_tr)$treatment)```

```pairwise.adonis(weighted_unifrac, phyloseq::sample_data(ps_tr)$treatment)```


_________________________________________________________________________________________________

## INSTALLATION

This repository is a fork of [pmartinezarbizu/pairwiseAdonis](https://github.com/pmartinezarbizu/pairwiseAdonis)
with bug fixes and compatibility updates for current versions of vegan (see NEWS.md in the package).
The instructions below install this fork. The package source lives in the `pairwiseAdonis/` subdirectory
of the repository, which is why the path has three parts.

The package is not on CRAN, so it is installed from GitHub with the `remotes` package
(`devtools::install_github()` works the same way). vegan, permute and cluster are installed
automatically if missing.

In your R session:

```
install.packages("remotes")
remotes::install_github("JWeld/pairwiseAdonis/pairwiseAdonis")
```

### Windows
Building from source requires Rtools, which can be installed from https://cran.r-project.org/bin/windows/Rtools/
(choose the version matching your R version). Then run the two lines above.

### Troubleshooting
If installation stops with "Error: Failed to install 'pairwiseAdonis' from GitHub: ... converted from warning"
(https://github.com/r-lib/remotes/issues/403), set this before installing:

```
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
```

To install the original upstream package instead, replace `JWeld` with `pmartinezarbizu` in the
`install_github()` call.

____________________________________
## Usage
```
library(pairwiseAdonis)
data(iris)
pairwise.adonis(iris[,1:4],iris$Species)

# For strata (blocks), following example of Jari Oksanen in adonis2. 
dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10,30)),field=gl(3,1) )
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(18)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(18)/2
Y <- data.frame(Agropyron, Schizachyrium)

pairwise.adonis2(Y ~ NO3/field, data = dat, strata = 'field')
```

for more examples see also
```?pairwise.adonis()```
```?pairwise.adonis2()```
_____________________________________________
## Citation

Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4.2
