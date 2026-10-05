#'@title Pairwise multilevel comparison using adonis2 accepting strata
#'
#'@description This is a wrapper function for multilevel pairwise comparison
#' using adonis2() from package 'vegan'. The function accepts interaction between factors and strata.
#'
#'@param x Model formula. The LHS is either community matrix or dissimilarity matrix (eg. from vegdist or dist;
#' a symmetric square matrix is also treated as a dissimilarity matrix, as in adonis2()).
#' See adonis2() for details. The RHS are factors that must be column names of a data.frame specified with argument data.
#' Pairwise comparisons are made between the levels of the first variable on the RHS, which should be a
#' factor, character or logical vector. A numeric variable gives a warning, because each pair of its
#' distinct values is compared (use e.g. \code{factor(x)} in the formula for numeric group codes).
#' If it is a factor, the pairs follow the order of its levels. The other terms of the formula are kept
#' in the model for each pair. \code{Condition()} terms for partial models (vegan >= 2.8-0) are passed
#' to adonis2() and are never used as grouping variable.
#'
#'@param data The data frame of independent variables having as column names the factors specified in formula.
#' It must have one row per row (observation) of the LHS.
#'
#'@param strata String. The name of the column with factors to be used as strata.
#'
#'@param nperm The number of permutations.
#'
#'@param p.adjust.m The p-value correction method, one of the methods supported by p.adjust(),
#' default is 'bonferroni'. The p-value in the first row of the table of each pair (see Details) is
#' adjusted for the number of pairwise comparisons and added as column \code{Pr(adj)}; the
#' significance codes then refer to the adjusted p-values. The unadjusted p-values stay in column
#' \code{Pr(>F)}. With 'none' no column is added.
#'
#'@param ... Any other parameter passed to adonis2, for example \code{by} (see Details) or
#' \code{na.action}.
#'
#'@details Since vegan 2.6-8 the default of adonis2() is \code{by = NULL}, an omnibus test of the whole
#' model for each pair. To test each term of a model with several terms
#' (e.g. \code{Y ~ NO3/field}) separately, as adonis2() did before vegan 2.6-8, pass
#' \code{by = "terms"} (or \code{by = "margin"}), which is forwarded to adonis2().
#'
#' The adjusted p-value refers to the first row of each table. This is the test of the grouping
#' variable when it is the only term of the model, or with \code{by = "terms"}, where the grouping
#' variable is the first term. With several terms and \code{by = NULL} the first row is the
#' omnibus test of the whole model for the pair.
#'
#' Observations with missing values in the variables on the RHS are dropped before the pairs are
#' formed when \code{na.action} (passed to adonis2) allows it, for example \code{na.action = na.omit};
#' with the default \code{na.fail} they give an error. Missing values in the strata column always
#' give an error.
#'
#'@return List of class "pwadstrata". The first element is the parent call; the remaining elements are the
#' anova tables returned by adonis2 for each unique pairwise combination of levels, with the
#' adjusted p-values in column \code{Pr(adj)} (unless \code{p.adjust.m = 'none'}).
#'
#'@author Pedro Martinez Arbizu
#'
#'@examples
#' data(iris)
#' pairwise.adonis2(iris[,1:4]~Species,data=iris)
#'
#' #For strata (blocks), Jari Oksanen recommends in the help of adonis2 to define the
#' #permutation matrix outside the adonis2 call
#' #In this example I have adapted the adonis2 example to have 3 factors in NO3
#'
#' dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10,30)),field=gl(3,1) )
#' Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(18)/2
#' Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(18)/2
#' Y <- data.frame(Agropyron, Schizachyrium)
#' perm <- how(nperm = 199)
#' setBlocks(perm) <- with(dat, field)
#' adonis2(Y ~ NO3, data = dat, permutations = perm)
#'
#'
#' # pairwise comparison
#' pairwise.adonis2(Y ~ NO3, data = dat, strata = 'field')
#' #notice the apostrophes in strata = 'field' !
#'
#' #this will give same results a doing adonis2 pairwise one by one
#' #(apart from the adjusted p-values in column Pr(adj))
#'
#' #for factors '0' and '10'
#' dat2 <- dat[dat$NO3 %in% c('0','10'),]
#' Y2 <- Y[dat$NO3 %in% c('0','10'),]
#' setBlocks(perm) <- with(dat2, field)
#' adonis2(Y2 ~ NO3, data = dat2, permutations = perm)
#' # and so on...
#'
#' # nested model, testing each term separately for each pair
#' pairwise.adonis2(Y ~ NO3/field, data = dat, strata = 'field', by = "terms")
#'
#'@export pairwise.adonis2
#'@importFrom utils combn
#'@importFrom vegan adonis2 vegdist
#'@import permute
#'@importFrom stats as.dist complete.cases model.frame na.fail na.pass p.adjust p.adjust.methods


pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, p.adjust.m = 'bonferroni', ... ) {

  if (!inherits(x, 'formula') || length(x) != 3L)
    stop("'x' must be a model formula with a response, e.g. Y ~ group")
  if (length(dim(data)) != 2L)
    stop("'data' must be a data frame")
  dots <- list(...)
  if ('permutations' %in% names(dots))
    stop("use arguments 'nperm' and 'strata' instead of 'permutations'")
  if (!is.null(strata) && !(is.character(strata) && length(strata) == 1 && strata %in% names(data)))
    stop("'strata' must be the name of a column of 'data', given as a string")
  p.adjust.m <- match.arg(p.adjust.m, p.adjust.methods)

#describe parent call function
ststri <- ifelse(is.null(strata),'Null',strata)
fostri <- as.character(x)
dotexpr <- match.call(expand.dots = FALSE)$...
dotstri <- ''
if (length(dotexpr)) {
  dotnames <- if (is.null(names(dotexpr))) rep('', length(dotexpr)) else names(dotexpr)
  dotstri <- paste0(', ', ifelse(dotnames == '', '', paste(dotnames, '= ')),
                    vapply(dotexpr, function(e) paste(deparse(e), collapse = ' '), ''),
                    collapse = '')
}

# environment of the formula: right hand side variables not found in 'data'
# and the response are looked up there, as in adonis2
  fenv <- environment(x)
  if (is.null(fenv)) fenv <- parent.frame()
#copy model formula
   x1 <- x
   environment(x1) <- fenv
# extract left hand side of formula
  lhs <- eval(x1[[2]], fenv, globalenv())
# a symmetric square matrix is a dissimilarity matrix (same rule as adonis2)
  if ((is.matrix(lhs) || is.data.frame(lhs)) && nrow(lhs) == ncol(lhs) &&
      is.numeric(as.matrix(lhs)) && isSymmetric(unname(as.matrix(lhs))))
    lhs <- as.dist(lhs)
  nobs <- if (inherits(lhs, 'dist')) attr(lhs, 'Size') else nrow(lhs)
  if (nobs != nrow(data))
    stop("the response and 'data' must have the same number of rows (observations): ",
         nobs, " vs ", nrow(data))

# right hand side only. Condition() terms are not candidates for the grouping
# variable, but their variables count for missing values
  x1[[2]] <- NULL
  rhs <- split_condition(x1[[2]])
  if (is.null(rhs$rest))
    stop("the formula needs a grouping variable outside Condition()")
  fgrp <- x1
  fgrp[[2]] <- rhs$rest
  fall <- x1
  fall[[2]] <- Reduce(function(a, b) call('+', a, b),
                      c(list(rhs$rest), lapply(rhs$cond, function(z) call('(', z))))
# create model.frame matrix (keep NA rows so that it stays aligned with data)
  rhs.frame <- model.frame(fall, data, na.action = na.pass)
  grp.frame <- if (length(rhs$cond)) model.frame(fgrp, data, na.action = na.pass) else rhs.frame
  if (!length(grp.frame))
    stop("the formula needs a grouping variable on the right hand side")

# missing values on the right hand side: the observations are dropped here, as
# adonis2 would do with its na.action, so that the response, data and strata
# stay aligned in every pair
  if (!all(complete.cases(rhs.frame))) {
    nafun <- if (is.null(dots$na.action)) na.fail else match.fun(dots$na.action)
    omit <- tryCatch(attr(nafun(rhs.frame), 'na.action'), error = function(e) e)
    if (inherits(omit, 'error'))
      stop("missing values in ",
           if (anyNA(grp.frame[[1]])) paste0("the grouping variable '", names(grp.frame)[1], "'")
           else "the right hand side variables",
           ": remove these observations or use na.action = na.omit")
    if (length(omit)) {
      keep <- !seq_len(nobs) %in% omit
      lhs <- if (inherits(lhs, 'dist')) as.dist(as.matrix(lhs)[keep, keep]) else lhs[keep, , drop = FALSE]
      data <- data[keep, , drop = FALSE]
      grp.frame <- grp.frame[keep, , drop = FALSE]
    }
  }
  if (!is.null(strata) && anyNA(data[[strata]]))
    stop("missing values in the strata column '", strata, "'")

# grouping variable: first variable on the right hand side
  gvar <- grp.frame[[1]]
  gname <- names(grp.frame)[1]
  if (!is.null(dim(gvar)))
    stop("the grouping variable '", gname, "' (the first variable on the right hand side) ",
         "must be a vector, not a matrix")
  grp <- as.character(gvar)
  lev <- if (is.character(gvar)) unique(grp) else levels(droplevels(as.factor(gvar)))
  if (length(lev) < 2)
    stop("the grouping variable '", gname, "' must have at least two levels")
  if (is.numeric(gvar))
    warning("the grouping variable '", gname, "' (the first variable on the right hand side) ",
            "is numeric: its ", length(lev), " distinct values are compared pairwise (",
            choose(length(lev), 2), " comparisons). Use factor(", gname,
            ") in the formula if these are group codes, or put the grouping factor first")

# create unique pairwise combination of factors
  co <- combn(lev, 2)

# create names vector
  nameres <- c('parent_call', paste(co[1, ], co[2, ], sep = '_vs_'))
#create results list
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres

#add parent call to res
res[['parent_call']] <- paste0(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ),
                               dotstri, ', p.adjust.m = ', p.adjust.m)

#start iteration trough pairwise combination of factors
 for(elem in seq_len(ncol(co))){

  idx <- grp %in% co[, elem]

#reduce model elements
	if(inherits(lhs,'dist')){
	    xred <- as.dist(as.matrix(lhs)[idx, idx])
	}else{
	    xred <- lhs[idx, , drop = FALSE]
	}

	mdat1 <- data[idx, , drop = FALSE]

# redefine formula: same right hand side, reduced response. adonis2 looks for
# the response in the environment of the formula (vegan <= 2.6) or in the
# calling frame (vegan >= 2.7), so 'xred' is kept in both: as a local variable
# here and in a child of the original formula environment, where functions and
# objects of the caller used on the right hand side are still found.
	xenv <- new.env(parent = fenv)
	xenv$xred <- xred
	xnew <- x
	xnew[[2]] <- as.name('xred')
	environment(xnew) <- xenv

#pass new formula to adonis2
	perm <- how(nperm = nperm)
	if(!is.null(strata)){
	    setBlocks(perm) <- mdat1[[strata]]
	    perm[['blocks.name']] <- strata
	}
	ad <- adonis2(xnew, data = mdat1, permutations = perm, ... )
# describe the pair instead of the internal call in the heading
	attr(ad, 'heading')[2] <- paste('Pairwise comparison', co[1, elem], 'vs', co[2, elem],
	                                'of', gname, 'in', paste(deparse(x), collapse = ' '))

  res[[nameres[elem+1]]] <- ad
  }

# adjust the p-values of the first row for the number of pairwise comparisons
  if (p.adjust.m != 'none') {
    p <- vapply(res[-1], function(a) a[1, 'Pr(>F)'], numeric(1))
    padj <- p.adjust(p, method = p.adjust.m)
    for (i in seq_along(padj)) {
      a <- res[[i + 1]]
      a[['Pr(adj)']] <- c(padj[i], rep(NA, nrow(a) - 1))
      attr(a, 'heading') <- c(attr(a, 'heading'),
                              paste0('Pr(adj): p-value of the first row adjusted for ', length(padj),
                                     ' pairwise comparisons (', p.adjust.m, ')\n'))
      res[[i + 1]] <- a
    }
  }
  class(res) <- c("pwadstrata", "list")
  return(res)
}

## split the right hand side of a formula into Condition() terms and the rest
split_condition <- function(e) {
  if (is.call(e) && identical(e[[1]], as.name('Condition')))
    return(list(rest = NULL, cond = list(e[[2]])))
  if (is.call(e) && identical(e[[1]], as.name('+')) && length(e) == 3L) {
    l <- split_condition(e[[2]])
    r <- split_condition(e[[3]])
    rest <- if (is.null(l$rest)) r$rest
            else if (is.null(r$rest)) l$rest
            else call('+', l$rest, r$rest)
    return(list(rest = rest, cond = c(l$cond, r$cond)))
  }
  list(rest = e, cond = list())
}


### Method print
#'@export
print.pwadstrata <- function(x, ...) {
  print(unclass(x), ...)
  invisible(x)
}

### Method summary
#'@export
summary.pwadstrata <- function(object, ...) {
  cat("Result of pairwise.adonis2:\n")
  cat("\n")
  print(unclass(object), ...)
  cat("\n")

  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(object)
}
