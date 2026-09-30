#'@title Pairwise multilevel comparison using adonis2
#'
#'@description This is a wrapper function for multilevel pairwise comparison
#' using adonis2() from package 'vegan'. The function returns adjusted p-values using p.adjust().
#'
#'@param x Data frame or matrix (the community table), or "dist" object (user-supplied distance matrix).
#' A symmetric square matrix is also treated as a distance matrix, as in adonis2().
#'
#'@param factors Vector (a column or vector with the levels to be compared pairwise).
#' Must have one entry per row (observation) of x.
#'
#'@param sim.function Function used to calculate the similarity matrix,
#' one of 'daisy' or 'vegdist' default is 'vegdist'. Ignored if x is a distance matrix.
#'
#'@param sim.method Similarity method from daisy or vegdist, default is 'bray'. Ignored if x is a distance matrix.
#'
#'@param p.adjust.m The p.value correction method, one of the methods supported by p.adjust(),
#' default is 'bonferroni'.
#'
#'@param reduce String. Restrict comparison to pairs including these factor levels. If more than one level,
#' separate by pipes like  reduce = 'setosa|versicolor'. Levels are matched exactly. Only the retained
#' comparisons are computed, and the p-value adjustment is applied to these comparisons only.
#'
#'@param perm The number of permutations, or a permutation design from permute::how().
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
#'
#'@export pairwise.adonis
#'@importFrom stats p.adjust as.dist
#'@importFrom utils combn
#'@importFrom vegan adonis2 vegdist
#'@importFrom cluster daisy


pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)
{
  factors <- as.character(factors)
  if (anyNA(factors))
    stop("'factors' contains missing values")

  ## a symmetric square matrix is a distance matrix (same rule as adonis2)
  if (!inherits(x, 'dist') && (is.matrix(x) || is.data.frame(x)) &&
      nrow(x) == ncol(x) && is.numeric(as.matrix(x)) &&
      isSymmetric(unname(as.matrix(x))))
    x <- as.dist(x)

  nobs <- if (inherits(x, 'dist')) attr(x, 'Size') else nrow(x)
  if (length(factors) != nobs)
    stop("'factors' must have one entry per row (observation) of 'x': ",
         length(factors), " vs ", nobs)

  lev <- unique(factors)
  if (length(lev) < 2)
    stop("'factors' must have at least two levels")

  co <- combn(lev, 2)

  ## restrict to pairs including the requested levels
  if (!is.null(reduce)) {
    keep.lev <- unlist(strsplit(reduce, '|', fixed = TRUE))
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
    if(inherits(x, 'dist')){
      x1 <- as.dist(as.matrix(x)[idx, idx])
    } else if (sim.function == 'daisy'){
      x1 <- daisy(x[idx, , drop = FALSE], metric = sim.method)
    } else {
      x1 <- vegdist(x[idx, , drop = FALSE], method = sim.method)
    }

    x2 <- data.frame(Fac = factors[idx])

    ad <- adonis2(x1 ~ Fac, data = x2, permutations = perm)
    pairs[elem] <- paste(co[1,elem],'vs',co[2,elem])
    Df[elem] <- ad$Df[1]
    SumsOfSqs[elem] <- ad$SumOfSqs[1]
    F.Model[elem] <- ad$F[1]
    R2[elem] <- ad$R2[1]
    p.value[elem] <- ad$`Pr(>F)`[1]
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  sig <- sig_codes(p.adjusted)

  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
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
