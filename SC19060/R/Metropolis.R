#' @Title Metropolis sampler
#' @name Metropolis
#' @description generate the Cauchy distribution by Metropolis sampler
#' @param sigma variance of the proposal distribution
#' @param x0 the initial value
#' @param N the number of samples
#' @return a randpm sample of Cauchy distribution and the number of candidate points rejected
#' @examples
#' \dontrun{
#' N<-2000
#' sigma<-.5
#' x0<-20
#' rw<-Metro.Cauchy(sigma,x0,N)
#' print(rw$k)
#' }
#' @export
Metro.Cauchy <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (1/(pi*(1+y^2))) / (1/(pi*(1+x[i-1]^2))))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      } 
  }
  return(list(x=x, k=k))
}