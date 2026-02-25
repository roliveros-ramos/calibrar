
#' Objective function between observed and simulated data
#'
#' Compute an objective value (error, negative log-likelihood, and/or penalty)
#' given observed data (`obs`) and model outputs (`sim`). This function is a thin
#' dispatcher: it resolves `FUN` via [match.fun()] and returns `FUN(obs = obs,
#' sim = sim, ...)`.
#'
#' The returned value is intended to be **minimised** by an optimiser: lower
#' values indicate a better match between `sim` and `obs` (or a lower penalty).
#'
#' @param obs Observed data as expected by `FUN`.
#'   Typically a numeric vector, matrix, or array. Missing values (`NA`) are
#'   permitted; see *Missing values* below.
#'
#' @param sim Simulated data matching `obs`, in the sense expected by `FUN`.
#'   For most pointwise criteria (e.g. `norm2`, `lnorm2`, `pois`), `sim` should
#'   have the same shape as `obs`. For composition criteria (e.g. `multinom`),
#'   `obs` and `sim` are expected to be matrices with rows representing samples
#'   and columns representing classes.
#'
#' @param FUN Objective function to apply. Can be:
#'   \itemize{
#'     \item a function, e.g. `FUN = norm2`;
#'     \item or a character string naming a function, e.g. `FUN = "norm2"`.
#'   }
#'   The function must accept arguments named `obs` and `sim` (additional
#'   arguments may be supplied via `...`) and must return a single numeric scalar.
#'
#' @param ... Additional arguments forwarded to `FUN`.
#'
#' @return A numeric scalar: the value of `FUN(obs = obs, sim = sim, ...)`.
#'   By convention this is an objective to be **minimised**.
#'
#' @section Conventions and expectations:
#' \itemize{
#'   \item **Minimisation**: all provided objectives are formulated so that lower
#'     values indicate a better fit (or weaker penalty).
#'   \item **Shapes**: `objFn()` does not reshape or recycle data. It is the
#'     caller's responsibility to supply `obs` and `sim` in a compatible form for
#'     the chosen `FUN`.
#'   \item **Scalar output**: `FUN` should return a length-one numeric value.
#'     Returning vectors is not supported by `objFn()` itself (although higher
#'     level workflows may aggregate vector-valued objectives elsewhere).
#' }
#'
#' @section Missing values:
#' Most built-in objectives stop if `all(is.na(obs))`. Otherwise, they ignore
#' missing values using `na.rm = TRUE` inside `sum()`/`mean()`. This means:
#' \itemize{
#'   \item partial `NA`s in `obs` are ignored;
#'   \item `NA`s in `sim` will also be dropped from sums when `na.rm = TRUE`
#'     (possibly masking simulation failures if not checked upstream).
#' }
#'
#' @section Numerical constraints and stability:
#' Some objectives impose additional constraints:
#' \itemize{
#'   \item **Poisson** (`pois`): uses `log(sim)`; `sim` must be strictly positive
#'     wherever `obs > 0` (otherwise `-Inf` may occur). Consider flooring the
#'     simulated intensity, e.g. `sim <- pmax(sim, 1e-12)`.
#'   \item **Log-scale** (`lnorm2`, `lnorm3`, `lnorm4`, `lnorm4b`): apply
#'     `log(x + tiny)`; both `obs + tiny` and `sim + tiny` must be positive.
#'     If your data may contain zeros, `tiny` should be chosen accordingly.
#'   \item **Compositions** (`multinom`): expects matrix inputs and uses row sums
#'     to convert counts/weights into proportions. See details below.
#' }
#'
#' @section Supported built-in objective functions:
#' The following functions are available with their current behaviour. All are
#' formulated as objectives to be minimised.
#'
#' \describe{
#'   \item{\code{norm2}}{
#'     Sum of squared errors on the original scale:
#'     \deqn{\sum (obs - sim)^2}
#'     Typical use: continuous observations with approximately additive errors.
#'   }
#'
#'   \item{\code{lnorm2}}{
#'     Sum of squared errors on the log scale:
#'     \deqn{\sum (\log(obs + tiny) - \log(sim + tiny))^2}
#'     Arguments: `tiny` (default `1e-2`), added before the log to avoid
#'     `log(0)`. Typical use: positive-valued data with multiplicative (lognormal)
#'     error structure.
#'   }
#'
#'   \item{\code{lnorm3}}{
#'     Log-scale squared error with an estimated multiplicative scaling factor
#'     \eqn{q}. Internally:
#'     \itemize{
#'       \item compute element-wise ratios `ratio <- obs/sim`;
#'       \item set `NaN` ratios to `NA`;
#'       \item estimate `q <- mean(ratio, na.rm = TRUE)`;
#'       \item return \eqn{\sum (\log(obs+tiny) - \log(sim+tiny) - \log(q))^2}.
#'     }
#'     Typical use: positive data where an overall multiplicative bias (scale
#'     mismatch) is expected and should not be fully penalised.
#'   }
#'
#'   \item{\code{lnorm4} / \code{lnorm4b}}{
#'     Extensions of `lnorm3` that add a penalty term to discourage extreme
#'     values of the scaling factor \eqn{q} (or extreme per-observation ratios).
#'     They rely on the helper \code{rangeq()}:
#'     \itemize{
#'       \item compute `ratio <- obs/sim`, estimate `q <- mean(ratio, na.rm=TRUE)`;
#'       \item compute a penalty using parameters `b` and `c`:
#'         \deqn{pen = n \cdot (\max(|\log_2(q)|, b)^c - b^c)}
#'         when `dump = TRUE` (used by `lnorm4`), or
#'         \deqn{pen = \sum (\max(|\log_2(ratio)|, b)^c - b^c)}
#'         when `dump = FALSE` (used by `lnorm4b`),
#'         where `n` is the number of non-missing ratios.
#'       \item add the penalty to the `lnorm3`-style objective.
#'     }
#'     Arguments: `tiny`, `b` (default `1`), `c` (default `2`).
#'     Typical use: log-scale fitting where scale drift must be controlled.
#'   }
#'
#'   \item{\code{pois}}{
#'     Poisson negative log-likelihood (up to constants):
#'     \deqn{-\sum (obs \log(sim) - sim)}
#'     Typical use: counts (or count-like rates) with Poisson observation error.
#'     Note: `sim` must be positive where `obs > 0`.
#'   }
#'
#'   \item{\code{multinom}}{
#'     A composition (multinomial-like) objective operating on matrices.
#'     Inputs are expected as `obs` and `sim` matrices with:
#'     \itemize{
#'       \item rows = samples (e.g. time steps, hauls, sites),
#'       \item columns = classes (e.g. age/size bins, categories).
#'     }
#'     Internal steps (high-level):
#'     \itemize{
#'       \item Let \eqn{A} be the number of classes (`A <- ncol(sim)`).
#'       \item Rows of `sim` that are all zeros (excluding rows that are all `NA`)
#'         are replaced by `1` on that row (interpreted as a uniform prior).
#'       \item Row sums are used to compute proportions:
#'         \eqn{Psim = sim/sum(sim)}, \eqn{Pobs = obs/sum(obs)} (row-wise).
#'       \item Rows with `sum(sim) == 0` are set to `NA` for numerical convenience.
#'       \item Rows with `sum(obs) == 0` are set to `NA` (interpreted as “no
#'         proportion data available”).
#'       \item A variance term `sigma2` and a small stabiliser `tiny` are used to
#'         define an objective that penalises discrepancies between `Pobs` and `Psim`.
#'     }
#'     Arguments: `size` (default `20`) and `tiny` (default `1e-3`).
#'     Interpretation: `size` plays the role of an effective sample size (larger
#'     values typically increase the weight of the composition fit).
#'   }
#'
#'   \item{\code{normp} / \code{re}}{
#'     Pure penalty on simulated values:
#'     \deqn{\sum sim^2}
#'     This ignores `obs` and can be used as a regulariser, or when `sim`
#'     represents a residual vector or deviates already computed upstream.
#'     `re` is an alias of `normp`.
#'   }
#'
#'   \item{\code{penalty}}{
#'     Scaled quadratic penalty:
#'     \deqn{n \cdot mean(sim^2)}
#'     Arguments: `n` (default `100`). This assumes a fixed sample size and can
#'     be used to put the penalty on a comparable scale across datasets.
#'   }
#' }
#'
#' @section Writing your own objective function:
#' You can supply any custom function via `FUN` provided it:
#' \itemize{
#'   \item accepts arguments named `obs` and `sim` (plus optional `...`);
#'   \item returns a length-one numeric scalar to be minimised;
#'   \item defines its own parameter checks and missing-value policy.
#' }
#'
#' @examples
#' ## Basic squared-error objective
#' obs <- c(1, 2, 3, NA, 5)
#' sim <- c(1.2, 1.9, 2.7, 4.0, 5.1)
#' objFn(obs, sim, FUN = "norm2")
#'
#' ## Log-scale objective (positive data)
#' obs <- c(0.1, 1, 10)
#' sim <- c(0.2, 0.9, 11)
#' objFn(obs, sim, FUN = lnorm2, tiny = 1e-2)
#'
#' ## Poisson objective (counts) with flooring for numerical safety
#' obs <- c(0, 3, 10, 2)
#' sim <- c(0, 2.5, 9.8, 1.9)
#' sim <- pmax(sim, 1e-12)
#' objFn(obs, sim, FUN = "pois")
#'
#' ## Composition objective (matrices: rows = samples, cols = classes)
#' obs <- rbind(c(10, 5, 0),
#'             c( 0, 0, 0),  # interpreted as “no composition data”
#'             c( 2, 1, 7))
#' sim <- rbind(c( 9, 6, 1),
#'             c( 0, 0, 0),  # replaced internally by sim+1 on that row
#'             c( 1, 2, 6))
#' objFn(obs, sim, FUN = "multinom", size = 20, tiny = 1e-3)
#'
#' ## Custom objective function
#' my_obj <- function(obs, sim, ...) {
#'   if (all(is.na(obs))) stop("All observed values are NA.")
#'   sum(abs(obs - sim), na.rm = TRUE)  # L1 error
#' }
#' objFn(obs = c(1, 2, NA), sim = c(1.1, 1.7, 3), FUN = my_obj)
#'
#' @seealso
#' \code{match.fun}
#'
#' @export
objFn = function(obs, sim, FUN, ...) {
  FUN = match.fun(FUN)
  output = FUN(obs=obs, sim=sim, ...)
  return(output)
}

#' @export
objFn = function(obs, sim, FUN, ...) {
  FUN = match.fun(FUN)
  output = FUN(obs=obs, sim=sim, ...)
  return(output)
}

#' @export
#' @rdname objFn 
fitness = objFn

# Penalties ---------------------------------------------------------------

normp = function(obs, sim, ...) {
  penalty = sum((sim)^2, na.rm=TRUE)
  return(penalty)
}

re = normp

penalty = function(obs, sim, n=100, ...) {
  # assumes a fixed sample size of 'n'
  penalty = n*mean((sim)^2, na.rm=TRUE)
  return(penalty)
}

# identical to penalty, used for testing.
penalty2 = penalty


# Likelihoods -------------------------------------------------------------

pois = function(obs, sim, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  nlogLike = -sum(obs*log(sim) - sim, na.rm=TRUE)
  return(nlogLike)
}

norm2 = function(obs, sim, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm2 = function(obs, sim, tiny=1e-2, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  obs = log(obs + tiny)
  sim = log(sim + tiny)
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm3  = function(obs, sim, tiny = 1e-2, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  ratio = obs/sim
  ratio[is.nan(ratio)] = NA
  q = mean(ratio, na.rm=TRUE)
  obs = log(obs+tiny) 
  sim = log(sim+tiny)
  nlogLike = sum((obs-sim-log(q))^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm4  = function(obs, sim, tiny = 1e-2, b=1, c=2, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  penq = rangeq(obs=obs, sim=sim, b=b, c=c, dump=TRUE)
  q = attr(penq, "q")
  obs = log(obs+tiny) 
  sim = log(sim+tiny)
  nlogLike = sum((obs-sim-log(q))^2, na.rm=TRUE)
  return(nlogLike + penq)
}

lnorm4b  = function(obs, sim, tiny = 1e-2, b=1, c=2, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  penq = rangeq(obs=obs, sim=sim, b=b, c=c, dump=FALSE)
  q = attr(penq, "q")
  obs = log(obs+tiny) 
  sim = log(sim+tiny)
  nlogLike = sum((obs-sim-log(q))^2, na.rm=TRUE)
  return(nlogLike + penq)
}

rangeq = function(obs, sim, b=1, c=2, dump=TRUE) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  ratio = obs/sim
  ratio[is.nan(ratio)] = NA
  n = sum(!is.na(ratio))
  q = mean(ratio, na.rm=TRUE)
  if(isTRUE(dump)) {
    pen = n*(pmax(abs((log2(q))), b)^c - b^c)
  } else {
    pen = sum((pmax(abs((log2(ratio))), b)^c - b^c))
  }
  attr(pen, "q") = q
  return(pen)
}

multinom = function(sim, obs, size=20, tiny=1e-3) {
  
  if(all(is.na(obs))) stop("All observed values are NA.")
  
  A = ncol(sim) # number of classes
  
  # checking for simulated values with only zeros and replace ones (assumed uniform)
  sim.allzero = apply(sim, 1, FUN=function(x) all(na.omit(x) == 0) & !all(is.na(x)) )
  sim[which(sim.allzero), ] = sim[which(sim.allzero), ] + 1
  
  sim.sum = rowSums(sim, na.rm=TRUE) # only zero for all NAs.
  obs.sum = rowSums(obs, na.rm=TRUE)
  
  # setting to NA for numerical convenience
  sim.sum[sim.sum==0] = NA
  
  # removing 'all zeros' from obs, because it means we have no proportion data
  obs.sum[which(obs.sum==0)] = NA  
  
  Psim     = sim/sim.sum
  Pobs     = obs/obs.sum
  
  sigma2 = ((1-Pobs)*Pobs + 1/A)/size
  
  error = log(exp(-((Pobs - Psim)^2)/(2*sigma2)) + tiny)
  
  nlogLike = -size*sum(error, na.rm=TRUE)
  return(nlogLike)
}
