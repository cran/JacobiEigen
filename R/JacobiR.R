##' The Jacobi Algorithm
##'
##' Eigenvalues and optionally, eigenvectore of a real symmetric matrix using the
##' classical Jacobi algorithm, (Jacobi, 1854)
##' @title The Jacobi Algorithm in Pure R
##' @param x a real symmetric matrix
##' @param only_values A logical value: Do you want eigenvalues only?
##' @param eps an error tolerance
##' @export JacobiR
##' @examples
##' (V <- crossprod(matrix(1:25, 5)))
##' JacobiR(V)
JacobiR <- function(x, only_values = FALSE,
                    eps = if(!only_values) .Machine$double.eps else
                      sqrt(.Machine$double.eps)) {
  n <- nrow(x)
  H <- if(only_values) NULL else diag(n)
  eps <- max(eps, .Machine$double.eps)

  if(n > 1) {
    lt <- which(lower.tri(x))

    repeat {
      k <- lt[which.max(abs(x[lt]))]  ## the matrix element
      j <- floor(1 + (k - 2)/(n + 1)) ## the column
      i <- k - n * (j - 1)            ## the row

      if(abs(x[i, j]) < eps) break

      Si <- x[, i]
      Sj <- x[, j]

      theta <- 0.5*atan2(2*Si[j], Sj[j] - Si[i])
      c <- cos(theta)
      s <- sin(theta)

      x[i, ] <- x[, i] <- c*Si - s*Sj
      x[j, ] <- x[, j] <- s*Si + c*Sj
      x[i,j] <- x[j,i] <- 0
      x[i,i] <- c^2*Si[i] - 2*s*c*Si[j] + s^2*Sj[j]
      x[j,j] <- s^2*Si[i] + 2*s*c*Si[j] + c^2*Sj[j]
      if(!only_values) {
        Hi <- H[, i]
        H[, i] <- c*Hi - s*H[, j]
        H[, j] <- s*Hi + c*H[, j]
      }
    }
  }
  list(values = as.vector(diag(x)), vectors = H)
}
