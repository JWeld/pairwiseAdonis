#'@title Pairwise multilevel comparison using adonis2 accepting strata
#'
#'@description This is a wrapper function for multilevel pairwise comparison
#' using adonis2() from package 'vegan'. The function accepts interaction between factors and strata.
#'
#'@param x Model formula. The LHS is either community matrix or dissimilarity matrix (eg. from vegdist or dist;
#' a symmetric square matrix is also treated as a dissimilarity matrix, as in adonis2()).
#' See adonis2() for details. The RHS are factors that must be column names of a data.frame specified with argument data.
#' Pairwise comparisons are made between the levels of the first variable on the RHS; the other
#' terms of the formula are kept in the model for each pair.
#'
#'@param data The data frame of independent variables having as column names the factors specified in formula.
#' It must have one row per row (observation) of the LHS.
#'
#'@param strata String. The name of the column with factors to be used as strata.
#'
#'@param nperm The number of permutations.
#'
#'@param ... Any other parameter passed to adonis2, for example \code{by} (see Details).
#'
#'@details Since vegan 2.6-8 the default of adonis2() is \code{by = NULL}, an omnibus test of the whole
#' model for each pair. To test each term of a model with several terms
#' (e.g. \code{Y ~ NO3/field}) separately, as adonis2() did before vegan 2.6-8, pass
#' \code{by = "terms"} (or \code{by = "margin"}), which is forwarded to adonis2().
#'
#'@return List of class "pwadstrata". The first element is the parent call; the remaining elements are the
#' anova tables returned by adonis2 for each unique pairwise combination of levels.
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
#'@importFrom stats as.dist model.frame na.pass


pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {

  if ('permutations' %in% names(list(...)))
    stop("use arguments 'nperm' and 'strata' instead of 'permutations'")
  if (!is.null(strata) && !(is.character(strata) && length(strata) == 1 && strata %in% names(data)))
    stop("'strata' must be the name of a column of 'data', given as a string")

#describe parent call function
ststri <- ifelse(is.null(strata),'Null',strata)
fostri <- as.character(x)

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
# create model.frame matrix (keep NA rows so that it stays aligned with data;
# adonis2 handles NA through its na.action argument)
  x1[[2]] <- NULL
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE, na.action = na.pass)

# grouping variable: first variable on the right hand side
  grp <- as.character(rhs.frame[,1])
  if (anyNA(grp))
    stop("missing values in the grouping variable '", names(rhs.frame)[1], "'")
  lev <- unique(grp)
  if (length(lev) < 2)
    stop("the grouping variable '", names(rhs.frame)[1], "' must have at least two levels")

# create unique pairwise combination of factors
  co <- combn(lev, 2)

# create names vector
  nameres <- c('parent_call', paste(co[1, ], co[2, ], sep = '_vs_'))
#create results list
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres

#add parent call to res
res[['parent_call']] <- paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm )

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

# redefine formula: same right hand side, reduced response. The response is
# put in a child of the original formula environment, so that functions and
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
	}
	ad <- adonis2(xnew, data = mdat1, permutations = perm, ... )

  res[[nameres[elem+1]]] <- ad
  }
  class(res) <- c("pwadstrata", "list")
  return(res)
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
