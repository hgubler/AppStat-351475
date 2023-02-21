#function to sample from a Gaussian mixture model
rbilnorm <- function(N, mu1, mu2, sigma1, sigma2, tau){
  ind <- I(runif(N) > tau)
  X <- rep(0,N)
  X[ind] <- rlnorm(sum(ind), mu1, sigma1)
  X[!ind] <- rlnorm(sum(!ind), mu2, sigma2)
  return(X)
}

#pdf from the Gaussian mixture model
dbilnorm <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dlnorm(x,mu1,sigma1) + tau*dlnorm(x,mu2,sigma2)
  return(y)
}

#EM algorithm for GMMM
EM = function(x, mu1, mu2, sigma1, sigma2, tau) {
  n = length(x)
  lold = 0
  lnew = 10
  while (abs(lold-lnew) > 10^-2) {
    #calculate log likelihood for convergence criteria
    lold = sum(log((1 - tau)*dlnorm(x, mean = mu1, sd = sqrt(sigma1)) + tau * dlnorm(x, mean = mu2, sd = sqrt(sigma2))))
    p = dlnorm(x, mean = mu2, sd = sqrt(sigma2)) * tau / dbilnorm(x, mu1, mu2, sqrt(sigma1), sqrt(sigma2), tau)
    #update the parameters according to ML estimator of slide 15, week 7
    tau = 1/n * sum(p)
    mu1 = sum(log(x) * (1 - p)) / sum(1 - p)
    mu2 = sum(log(x) * p) / sum(p)
    sigma1 = sum((1 - p) * ((log(x) - mu1)^2)) / sum(1 - p)
    sigma2 = sum(p * ((log(x) - mu2)^2)) / sum(p)
    #calculate log likelihood for convergence criteria
    lnew = sum(log((1 - tau)*dlnorm(x, mean = mu1, sd = sqrt(sigma1)) + tau * dlnorm(x, mean = mu2, sd = sqrt(sigma2))))
  }
  return(c(mu1, mu2, sigma1, sigma2, tau))
}

data = read.csv(file = "Project-1/1_snow_particles.csv")
#Check histogram of the data
hist(data$retained....)
#jittering the data
n = length(data$X)
noise = runif(n, min = 0, max = 0.05)
jittered_data = data$startpoint + noise
hist(jittered_data, breaks = 20)
testdata = rbilnorm(10000, 1, 5, 1, 1, 0.6)
#define starting parameters of the bilognormal model
mu1 <- 1
mu2 <- 3
sigma1 <- 0.5
sigma2 <- 1
tau <- 0.5
estim = EM(testdata, mu1, mu2, sigma1, sigma2, tau)
print(estim)

