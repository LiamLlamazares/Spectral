#' Numerically stable sinc function, \code{sin(x)/x}
#'
#' @param x A numeric vector
#'
#' @return `sin(x)/x`
#' @export
#'
#' @examples
sinc <- function(x) {
  f <- sin(x)
  # Use Taylor expansion close to x=0
  # Minimiser of the relative error:
  # sin(x)/x = 1-x^2/6+x^4/120-O(x^6)...
  # Truncated Taylor: 1-x^2/6
  # Approximate relative computational error for sin(x)/x :
  #   |(1 + eps_1)/(1 + eps_2) - 1| <= 2 * eps
  # Set errors equal (approximatively):
  # x^4 / 120 = 2 * eps
  # x = (240 * eps)^{1/4}
  x0 <- (240 * .Machine$double.eps)^(1/4)
  ok <- abs(x) > x0
  if (all(ok)) {
    f <- f / x
  } else {
    f[ok] <- f[ok] / x[ok]
    f[!ok] <- 1 - x[!ok]^2 / 6
  }
  f
}

#' Folded spectrum
#'
#' @param omega
#' @param S_fun
#' @param h
#' @param pointwise
#' @param n_folds
#' @param fold If `TRUE`, computes a folded spectrum Stilde. That is, the sum of
#' If `FALSE`, computes
#' a truncated spectrum.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fold_spectrum <- function(omega, S_fun, h, pointwise = TRUE, n_folds = 10,
                          fold = TRUE,
                          ...) {
  S <- 0
  done <- FALSE
  if (fold) {
    kk <- seq(-n_folds, n_folds, by = 1)
  } else {
    kk <- 0
  }
  k_ <- rep(1, length(h))
  while (!done) {
    # Shifts the frequencies
    omega_ <- omega + kronecker(
      as.matrix(rep(1, NROW(omega))),
      t(2 * pi * kk[k_] / h)
    )
    # Scalin g by default is 1
    scaling <- 1
    if (!pointwise) {
      for (d in seq_along(h)) {
        scaling <- scaling * sinc(omega_[, d] * h[d] / 2)^2
      }
    }
    # Integration step
    S <- S + S_fun(omega_, ...) * scaling

    #Updates for next iteration
    idx <- which(k_ < length(kk))
    if (length(idx) == 0) {
      done <- TRUE
    } else {
      idx <- max(idx)
      k_[idx] <- k_[idx] + 1
      if (idx < length(h)) {
        k_[(idx + 1):length(h)] <- 1
      }
    }
  }
  global_scaling <- 1
  for (d in seq_along(h)) {
    global_scaling <- global_scaling *
      (omega[, d] >= -pi / h[d]) *
      (omega[, d] < pi / h[d])
  }
  S <- S * global_scaling
  S
}

#' Shift needed to be able to apply the Fast Fourier tranform
#'
#' @param x
#' @param dim integer vector of dimension extents, suitable for `array()`.
#' The length of `dim` is the dimension for the intended fft transformation
#'
#' @return
#' @export
#'
#' @examples
fftshift <- function(x, dim = length(x)) {
  stopifnot(prod(dim) == length(x))
  x <- array(as.vector(x), dim = dim)
  if (length(dim) == 1) {
    x <- x[
      c(seq_len(dim[1] / 2) + dim[1] / 2, seq_len(dim[1] / 2)),
      drop = FALSE
    ]
  } else if (length(dim) == 2) {
    x <- x[
      c(seq_len(dim[1] / 2) + dim[1] / 2, seq_len(dim[1] / 2)),
      c(seq_len(dim[2] / 2) + dim[2] / 2, seq_len(dim[2] / 2)),
      drop = FALSE
    ]
  } else {
    stop(paste0("Dimension ", length(dim), " not implemented."))
  }
  x
}
#' Calculates the covariance function from the spectrum
#'
#' @param S
#' @param dim
#' @param h
#'
#' @return
#' @export
#'
#' @examples
S2C <- function(S, dim, h) {
  if (length(dim) == 1) {
    fft <- fftwtools::fftw
  } else if (length(dim) == 2) {
    fft <- fftwtools::fftw2d
  } else {
    stop(paste0("Dimension ", length(dim), " is not implemented."))
  }
  C <- fftshift(Re(fft(fftshift(S, dim))), dim)
  C * prod(2 * pi / h / dim)
}

#' Spectral sampling
#'
#' @param S Evaluated spectrum
#' @param dim Number of frequency space integration points for each dimension
#' @param h Spatial step sizes
#' @param seed The random seed
#' @param conjugate If `TRUE`, use the complex conjugate property of the
#'   spectral process to generate a single realisation. If `FALSE`, the
#'   real and imaginary parts will be two independent realisations
#'
#' @return
#' @export
#'
#' @examples

S2sample <- function(S, dim, h, seed = NULL, conjugate = FALSE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (length(dim) == 1) {
    fft <- fftwtools::fftw
  } else if (length(dim) == 2) {
    fft <- fftwtools::fftw2d
  } else {
    stop(paste0("Dimension ", length(dim), " is not implemented."))
  }
  SD <- sqrt(S * prod(2 * pi / h / dim))
  z <- (rnorm(prod(dim)) + 1i*rnorm(prod(dim))) * SD
  if (conjugate) {
    k <- expand.grid(lapply(dim, function(d) seq_len(d) - 1 - d/2))
    # Find symmetry pairs and make them complex conjugate
    if (length(dim) == 1) {
      # Note: this can be implemented much more efficiently...
      pair <- rep(NA_integer_, nrow(k))
      for (idx in seq_len(nrow(k))) {
        pair_ <- k[idx,1] == -k[,1]
        if (any(pair_)) {
          pair[idx] <- which(pair_)
        }
      }
      idx_ <- (k[,1] < 0) & !is.na(pair)
      z[idx_] <- Re(z[pair[idx_]]) - 1i * Im(z[pair[idx_]])
      idx_ <- (k[,1] == 0) | is.na(pair)
      z[idx_] <- Re(z[idx_])
    } else { # length(dim) == 2
      # Note: this can be implemented much more efficiently...
      pair <- rep(NA_integer_, nrow(k))
      for (idx in seq_len(nrow(k))) {
        pair_ <- (k[idx,1] == -k[,1]) & (k[idx,2] == -k[,2])
        if (any(pair_)) {
          pair[idx] <- which(pair_)
        }
      }
      idx_ <- (k[,1] < 0) & !is.na(pair)
      z[idx_] <- Re(z[pair[idx_]]) - 1i * Im(z[pair[idx_]])
      idx_ <- (k[,1] == 0) & (k[,2] < 0) & !is.na(pair)
      z[idx_] <- Re(z[pair[idx_]]) - 1i * Im(z[pair[idx_]])
      idx_ <- ((k[,1] == 0) & (k[,2] == 0)) | is.na(pair)
      z[idx_] <- Re(z[idx_])
    }
  } else { # !conjugate
    k <- expand.grid(lapply(dim, function(d) seq_len(d) - 1 - d/2))
    # Find lower edges and set to zero
    if (length(dim) == 1) {
      edge <- k[,1] == -dim[1]/2
      z[edge] <- 0
    } else { # length(dim) == 2
      edge <- (k[,1] == -dim[1]/2) | (k[,2] == -dim[2]/2)
      z[edge] <- 0
    }
  }

  z <- fftshift(z, dim)
  if (length(dim) == 1) {
    sample <- (fft(z, dim))
  } else {
    sample <- (fft(matrix(z, dim[1], dim[2])))
  }
  sample
}



#' Title
#'
#' @param dim
#' @param L
#'
#' @return
#' @export
#'
#' @examples
make_x <- function(dim, L) {
  x <- list()
  for (d in seq_along(dim)) {
    x[[paste0("x", d)]] <-
      (seq_len(dim[d]) - 1 - dim[d] / 2) / dim[d] * L[d]
  }
  x
}

#' Locations for sampling results
#'
#' @param dim Number of spectral integration points, for each dimension
#' @param h The spatial step lengths, for each dimension
#'
#' @return
#' @export
#'
#' @examples
make_x_sampling <- function(dim, h) {
  L <- dim * h
  x <- list()
  for (d in seq_along(dim)) {
    x[[paste0("x", d)]] <-
      (seq_len(dim[d]) - 1) / dim[d] * L[d]
  }
  x
}

#' Frequencies for covariance calculations
#'
#' @param dim
#' @param L
#'
#' @return
#' @export
#'
#' @examples
make_omega <- function(dim, L) {
  h <- L / dim
  w <- list()
  for (d in seq_along(dim)) {
    w[[paste0("w", d)]] <-
      seq(-(dim[d] / 2), dim[d] / 2 - 1, by = 1) / (dim[d] / 2) * pi / h[d]
  }
  w
}

#' Frequencies for sampling
#'
#' @param dim Number of spectral integration points, for each dimension
#' @param h The spatial step lengths, for each dimension
#'
#' @return
#' @export
#'
#' @examples
make_omega_sampling <- function(dim, h) {
  w <- list()
  for (d in seq_along(dim)) {
    w[[paste0("w", d)]] <-
      seq(-(dim[d] / 2), dim[d] / 2 - 1, by = 1) / (dim[d] / 2) * pi / h[d]
  }
  w
}

## Code specific to anisotropic implementation

#' H matrix
#'
#' @param v A two by two constant vector

#' @return A 2x2 matrix H
#' @export
#'
#' @examples
#' v <- c(-1,0)
#' H <-calcH(v)

calcH <- function(v){
  # Pre-calculate common terms
  v_magnitude <- sqrt(v[1]^2+v[2]^2)
  cosh_val <- cosh(v_magnitude)
  sinh_over_magnitude <- sinh(v_magnitude) / v_magnitude

  # Calculate H
  H <- matrix(0, nrow = 2, ncol = 2)
  H[1, ] <- c(cosh_val + sinh_over_magnitude * v[1], sinh_over_magnitude * v[2])
  H[2, ] <- c(sinh_over_magnitude * v[2], cosh_val - sinh_over_magnitude * v[1])
  return(as.matrix(H))
}

#' Spectrum of anisotropic field
#'
#' @param v A two by two constant vector
#' @param logkappa The logarithm of the inverse correlation range
#' @param omega Frequency in R^2 at which spectrum is calculated

#' @return The spectrum S(omega) of the anisotropic field with parameters (kappa,H,sigma )
#' @export
#'
#' @examples
#' logkappa <- 2
#' v <- c(1,2)
#' omega <- c (1,1)
#' #' S_aniso(v = v, logkappa = logkappa, omega = omega)

S_aniso <- function(logkappa, v,sigma, omega) {
  #Calcuates kappa and H
  kappa <- exp(logkappa)
  H <-calcH(v)

  # Calculate omega^T H omega
  omega_H_omega <- as.numeric(t(omega) %*% H %*% omega)

  # Calculate the spectrum
  result <- 4*pi*kappa*sigma^2 / ((2 * pi)^2) / (kappa^2 + omega_H_omega)^2

  return(result)
}
#' Spectrum of anisotropic field applied to matrix of omegas
#'
#' @param v A two by two constant vector
#' @param logkappa The logarithm of the inverse correlation range
#' @param omega Matrix of size mx2 (w1,...,wm)

#' @return A vector of size n (S(w1),...,S(wm))
#' Where S is the spectrum of the anisotropic field with parameters (kappa,H )
#' @export
#'
#' @examples
#' logkappa <- 2
#' v <- c(1,2)
#' n <- c(256, 256)
#' L <- c(1, 1)
#' omega_ <- make_omega(n, L)
#' omega <- as.matrix(expand.grid(omega_))
#' S_aniso_matrix(logkappa=logkappa, v=v, omega=omega)



S_aniso_matrix <- function(logkappa, v, omega, sigma =1) {
  # Apply S_aniso to each row of omega_matrix
  new <-function(omega) S_aniso(logkappa=logkappa, v=v, sigma, omega=omega)
  apply(omega, 1, new)
}





