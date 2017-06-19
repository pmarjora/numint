

splitN <- function(N, nchunks) {
  nchunks <- min(N, nchunks)
  ans <- rep(N %/% nchunks, nchunks)
  ans[nchunks] <- ans[nchunks] + N %% nchunks
  ans
}

#' Monte Carlo integration
#' This function performes numerical integration using monte carlo integration
#' @param fn is a function
#' @param ... are further arguments to be passed to fn
#' @param a is the lower bound
#' @param b is the upper bound
#' @param N is the sample size
#' @param ncores is the number of processors to call
#' @param cl is an object of class cluster
#' @export
#'
#' @family Numerical methods
#' @references Boyle, E. A., Li, Y. I., & Pritchard, J. K. (2017). An Expanded View of Complex Traits: From Polygenic to Omnigenic. Cell, 169(7), 1177???1186. http://doi.org/10.1016/j.cell.2017.05.038
#'
#'
#' @return An object of class "numint"
#' @aliases NumInt
num_int <- function(fn, ..., a, b, N = 100, ncores = 1, cl = NULL) {

  # Getting the call
  call <- match.call()

  # Checking length
  if (length(b) != length(a))
    stop("-a- and -b- must have the same length.")

  # Checking values
  if (!all(a < b))
    stop("There are some values a > b (all must a < b)")

  # Checking parallel
  if (!length(cl)) {
    cl <- makeCluster(ncores)
    toload <- loadedNamespaces()
    invisible(parallel::clusterCall(cl, function(x) {
      sapply(x, library, character.only = TRUE)
    }, x = toload))
    on.exit(parallel::stopCluster(cl))
  }

  # Sampling
  # samp <- Map(function(lb, ub) stats::runif(N, lb, ub), lb = a, ub = b)
  # samp <- do.call(cbind, samp)

  # Computing the volume
  V <- prod(b - a)

  # Computing density
  f       <- function(x) fn(x, ...)

  # Distributing indices across processors
  # idx     <- splitIndices(N, length(cl))
  # fsample <- parLapply(cl, lapply(idx, function(w) samp[w,,drop=FALSE]), f)
  fsample <- parallel::parLapply(cl, splitN(N, length(cl)), function(n, a, b) {
    samp <- Map(function(lb, ub) runif(n, lb, ub), lb = a, ub = b)
    samp <- do.call(cbind, samp)
    f(samp)
  }, a = a, b = b)
  fsample <- unlist(fsample)

  ans     <- V*sum(fsample)/N

  # Preparing arguments
  env <- new.env()
  environment(f) <- env
  args <- list(...)
  environment(args) <- env

  # Returning an object of class numint
  structure(
    list(
      val = ans,
      vol = V,
      fsample = fsample,
      call = call,
      sd = sd(fsample*V),
      N=N,
      a = a,
      b = b,
      f = f,
      args = args
    ),
    class = "numint"
  )
}

#'  Plotting method for num_int
#' @param x the x-coordinate
#' @param y the y-coordinate (ignored, but need tbecause the standard plot function has them)
#' @param main the title
#' @param col the color to use
#' @param ...  additional parameters to pass
#' @return it returns the coordinates
#' @export
plot.numint <- function(x, y = NULL, main = "Monte Carlo Integration", col=blues9[4],...) {

  n <- 100

  # Computing values
  xran <- c(x$a[1], x$b[1])
  vals <- NULL
  if (length(x$a > 1)) {
    vals <- Map(function(a,b) rep((a + b)/2, n), a = x$a[-1], x$b[-1])
    vals <- do.call(cbind, vals)
  }

  # Computing coordinates
  vals <- cbind(seq(xran[1], xran[2], length.out = n), vals)
  y <- apply(vals, 1, x$f)


  # Adding missing points
  coordinates <- cbind(vals[,1], y)
  coordinates <- rbind(coordinates, c(xran[2], 0), c(xran[1], 0))

  # Plotting
  plot.new()
  plot.window(xlim = xran, ylim = range(y))
  polygon(coordinates, col = col, ...)

  # Adding axis
  axis(1);axis(2)

  # A nice title
  title(main = main)

  # And a neat legend
  legend("topright",
         legend = substitute(
           Volume~~a %+-% b,
           list(
             a = sprintf("%.4f", x$val),
             b = sprintf("%.4f", x$sd))
         ),
         bty = "n"
  )

  # Returning the coordinates used for the plot
  invisible(cbind(x = vals, y = y))
}

# printing method
print.numint <- function(x, ...) {
  with(x, cat(sprintf("MONTE CARLO INTEGRATION\nN: %i\nVolume: %.4f\n", N, vol)))
  with(x, cat(sprintf("%.4f +- %.4f", val, sd)))

  invisible(x)
}


#' Addnums
#' @rdname addnums
#' @export
#' @param x An object of class \code{funnypkg_addnums}.
#' @param y Ignored.
#' @param ... Further arguments passed to
#' \code{\link[graphics:plot.window]{plot.window}}.
plot.funnypkg_addnums <- function(x, y = NULL, ...) {
  graphics::plot.new()
  graphics::plot.window(xlim = range(unlist(c(0,x))), ylim = c(-.5,1))
  graphics::axis(1)
  with(x, graphics::segments(0, 1, ab, col = "blue", lwd=3))
  with(x, graphics::segments(0, 0, a, col = "green", lwd=3))
  with(x, graphics::segments(a, .5, a + b, col = "red", lwd=3))
  graphics::legend("bottom", col = c("blue", "green", "red"),
                   legend = c("a+b", "a", "b"), bty = "n",
                   ncol = 3, lty = 1, lwd=3)
}
