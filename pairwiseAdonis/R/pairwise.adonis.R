#'@title Pairwise multilevel comparison using adonis2
#'
#'@description This is a wrapper function for multilevel pairwise comparison
#' using adonis2() from package 'vegan'. The function returns adjusted p-values using p.adjust().
#'
#'@param x Data frame or matrix (the community table), or "dist" object (user-supplied distance matrix).
#' A symmetric square matrix is also treated as a distance matrix, as in adonis2().
#'
#'@param factors Vector (a column or vector with the levels to be compared pairwise).
#' Must have one entry per row (observation) of x. If it is a factor, the pairs follow the order
#' of its levels, otherwise the order in which the values first appear.
#'
#'@param sim.function Function used to calculate the similarity matrix,
#' one of 'daisy' or 'vegdist' default is 'vegdist'. Ignored if x is a distance matrix.
#'
#'@param sim.method Similarity method from daisy or vegdist, default is 'bray'. With daisy it must be
#' one of 'euclidean', 'manhattan' or 'gower'. Ignored if x is a distance matrix.
#'
#'@param p.adjust.m The p.value correction method, one of the methods supported by p.adjust(),
#' default is 'bonferroni'.
#'
#'@param reduce String. Restrict comparison to pairs including these factor levels. If more than one level,
#' separate by pipes like  reduce = 'setosa|versicolor'. Levels are matched exactly (spaces around the
#' pipes are ignored). Only the retained comparisons are computed, and the p-value adjustment is
#' applied to these comparisons only.
#'
#'@param perm The number of permutations, or a permutation design from permute::how().
#' Blocks defined in the design (setBlocks) are reduced to the observations of each pair.
#' Designs with plot-level strata (setPlots), restricted permutations within blocks
#' (Within() of type 'series' or 'grid') and permutation matrices are not accepted,
#' because each pair uses a different subset of the observations.
#'
#'@param ... Further named arguments. \code{sqrt.dist}, \code{add} and \code{parallel} are passed
#' to adonis2(), all others to the distance function (vegdist() or daisy()), for example
#' \code{binary = TRUE} for presence/absence indices with vegdist(). Arguments for the distance
#' function are ignored if x is a distance matrix.
#'
#'@return Table (data frame of class "pwadonis") with the pairwise factors, Df, SumsOfSqs, F-values, R^2,
#' p.value, adjusted p.value and significance codes
#' (0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1, applied to the adjusted p-values).
#'
#'@author Pedro Martinez Arbizu & Sylvain Monteux
#'
#'@examples
#' data(iris)
#' pairwise.adonis(iris[,1:4],iris$Species)
#'
#'
#'#similarity euclidean from vegdist and holm correction
#' pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='vegdist',
#' sim.method='euclidian',p.adjust.m='holm')
#'
#'#identical example using a distance matrix as an input
#' dist_matrix=vegdist(iris[,1:4],method="euclidean")
#' pairwise.adonis(dist_matrix,factors=iris$Species,
#' p.adjust.m='holm')
#'
#'#similarity manhattan from daisy and bonferroni correction
#' pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='daisy',
#' sim.method='manhattan',p.adjust.m='bonferroni')
#'
#'#Restrict comparison to only some factors
#'pairwise.adonis(iris[,1:4],iris$Species, reduce='setosa')
#'
#'#for more than one factor separate by pipes
#'pairwise.adonis(iris[,1:4],iris$Species, reduce='setosa|versicolor')
#'
#'#presence/absence (Jaccard) dissimilarities: binary = TRUE is passed to vegdist()
#' data(dune, dune.env)
#' pairwise.adonis(dune, dune.env$Management, sim.method = 'jaccard', binary = TRUE)
#'
#'
#'@export pairwise.adonis
#'@importFrom stats p.adjust p.adjust.methods as.dist
#'@importFrom utils combn
#'@importFrom vegan adonis2 vegdist
#'@importFrom cluster daisy
#'@importFrom permute getBlocks getStrata getType


pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999, ...)
{
  if (anyNA(factors))
    stop("'factors' contains missing values")
  ## pairs follow the order of the levels of a factor, else the order of appearance
  lev <- if (is.factor(factors)) levels(droplevels(factors)) else unique(as.character(factors))
  factors <- as.character(factors)
  sim.function <- match.arg(sim.function, c('vegdist', 'daisy'))
  p.adjust.m <- match.arg(p.adjust.m, p.adjust.methods)

  ## further arguments: some for adonis2, the others for the distance function
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(names(dots) == '')))
    stop("arguments in '...' must be named")
  if (any(c('permutations', 'strata') %in% names(dots)))
    stop("use argument 'perm' for the number of permutations or a permutation design")
  ad.args <- dots[names(dots) %in% c('sqrt.dist', 'add', 'parallel')]
  dist.args <- dots[!names(dots) %in% names(ad.args)]

  ## a symmetric square matrix is a distance matrix (same rule as adonis2)
  if (!inherits(x, 'dist') && (is.matrix(x) || is.data.frame(x)) &&
      nrow(x) == ncol(x) && is.numeric(as.matrix(x)) &&
      isSymmetric(unname(as.matrix(x))))
    x <- as.dist(x)

  if (!inherits(x, 'dist') && sim.function == 'daisy' &&
      !sim.method %in% c('euclidean', 'manhattan', 'gower'))
    stop("with sim.function = 'daisy', sim.method must be one of ",
         "'euclidean', 'manhattan' or 'gower', not '", sim.method, "'")

  nobs <- if (inherits(x, 'dist')) attr(x, 'Size') else nrow(x)
  if (length(factors) != nobs)
    stop("'factors' must have one entry per row (observation) of 'x': ",
         length(factors), " vs ", nobs)

  if (length(lev) < 2)
    stop("'factors' must have at least two levels")

  ## permutation design: blocks are reduced to each pair, other
  ## observation-specific designs cannot be reduced
  blocks <- NULL
  if (is.matrix(perm))
    stop("'perm' cannot be a permutation matrix: each pair uses a subset of the observations")
  if (inherits(perm, 'how')) {
    if (!is.null(getStrata(perm, which = 'plots')))
      stop("permutation designs with plot-level strata are not supported: ",
           "use blocks (permute::setBlocks) or pairwise.adonis2() with 'strata'")
    within.type <- getType(perm, which = 'within')
    if (!within.type %in% c('free', 'none'))
      stop("restricted permutations within blocks (type '", within.type, "') are not supported: ",
           "their spatial or temporal structure does not hold for the subset of observations in each pair")
    blocks <- getBlocks(perm)
    if (!is.null(blocks) && length(blocks) != nobs)
      stop("blocks of 'perm' must have one entry per observation of 'x'")
    if (anyNA(blocks))
      stop("blocks of 'perm' contain missing values")
  }

  co <- combn(lev, 2)

  ## restrict to pairs including the requested levels
  if (!is.null(reduce)) {
    keep.lev <- trimws(unlist(strsplit(reduce, '|', fixed = TRUE)))
    unknown <- setdiff(keep.lev, lev)
    if (length(unknown))
      warning("level(s) in 'reduce' not found in 'factors': ",
              paste(unknown, collapse = ', '))
    co <- co[, co[1, ] %in% keep.lev | co[2, ] %in% keep.lev, drop = FALSE]
    if (ncol(co) == 0)
      stop("no pairwise comparison left after applying 'reduce'")
  }

  pairs <- character(ncol(co))
  Df <- SumsOfSqs <- F.Model <- R2 <- p.value <- numeric(ncol(co))

  for(elem in seq_len(ncol(co))){
    idx <- factors %in% co[, elem]
    ## calls are built with do.call() on names, so that the data are not
    ## copied into the calls stored by vegdist() and adonis2()
    if(inherits(x, 'dist')){
      x1 <- as.dist(as.matrix(x)[idx, idx])
    } else {
      xsub <- x[idx, , drop = FALSE]
      x1 <- if (sim.function == 'daisy')
        do.call('daisy', c(list(quote(xsub), metric = sim.method), dist.args))
      else
        do.call('vegdist', c(list(quote(xsub), method = sim.method), dist.args))
    }

    x2 <- data.frame(Fac = factors[idx])

    perm.i <- perm
    if (!is.null(blocks))
      setBlocks(perm.i) <- blocks[idx]

    ad <- do.call('adonis2', c(list(quote(x1 ~ Fac), data = quote(x2),
                                    permutations = quote(perm.i)), ad.args))
    pairs[elem] <- paste(co[1,elem],'vs',co[2,elem])
    Df[elem] <- ad$Df[1]
    SumsOfSqs[elem] <- ad$SumOfSqs[1]
    F.Model[elem] <- ad$F[1]
    R2[elem] <- ad$R2[1]
    p.value[elem] <- ad$`Pr(>F)`[1]
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  sig <- sig_codes(p.adjusted)

  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig,
                          stringsAsFactors = FALSE)
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}

## standard significance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
sig_codes <- function(p) {
  sig <- rep('', length(p))
  sig[!is.na(p) & p <= 0.1] <- '.'
  sig[!is.na(p) & p <= 0.05] <- '*'
  sig[!is.na(p) & p <= 0.01] <- '**'
  sig[!is.na(p) & p <= 0.001] <- '***'
  sig
}

### Method summary
#'@export
summary.pwadonis <- function(object, ...) {
  cat("Result of pairwise.adonis:\n")
  cat("\n")
  print.data.frame(object, ...)
  cat("\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(object)
}
## end of method summary
