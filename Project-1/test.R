#function to sample from a Gaussian mixture model
rbilnorm <- function(N, mu1, mu2, sigma1, sigma2, tau){
  ind <- I(runif(N) > tau)
  X <- rep(0,N)
  X[ind] <- rlnorm(sum(ind), mu1, sigma1)
  X[!ind] <- rlnorm(sum(!ind), mu2, sigma2)
  return(X)
}

#pdf from the Gaussian mixture model
dbilnorm2 <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dlnorm(x,mu1,sigma1) + tau*dlnorm(x,mu2,sigma2)
  return(y)
}

dlnorm_test = function(x, mu, sigma) {
  y = 1 / (x * sigma * sqrt(2 * pi)) * exp(-(log(x) - mu)^2 / (2 * sigma^2))
  return(y)
}
dbilnorm <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dlnorm_test(x,mu1,sigma1) + tau*dlnorm_test(x,mu2,sigma2)
  return(y)
}

#EM algorithm for GMMM
EM = function(x, mu1, mu2, sigma1, sigma2, tau) {
  n = length(x)
  lold = 0
  lnew = 10
  while (abs(lold-lnew) > 10^-2) {
    #calculate log likelihood for convergence criteria
    lold = sum(log((1 - tau)*dlnorm(x, mean = mu1, sd = sigma1) + tau * dlnorm(x, mean = mu2, sd = sigma2)))
    p = dlnorm(x, mean = mu2, sd = sigma2) * tau / dbilnorm(x, mu1, mu2, sigma1, sigma2, tau)
    #update the parameters according to ML estimator of slide 15, week 7
    tau = 1/n * sum(p)
    mu1 = sum(log(x) * (1 - p)) / sum(1 - p)
    mu2 = sum(log(x) * p) / sum(p)
    sigma1 = sqrt(sum((1 - p) * ((log(x) - mu1)^2)) / sum(1 - p))
    sigma2 = sqrt(sum(p * ((log(x) - mu2)^2)) / sum(p))
    #calculate log likelihood for convergence criteria
    lnew = sum(log((1 - tau)*dlnorm(x, mean = mu1, sd = sigma1) + tau * dlnorm(x, mean = mu2, sd = sigma2)))
  }
  return(c(mu1, mu2, sigma1, sigma2, tau))
}

data = read.csv(file = "Project-1/1_snow_particles.csv")
num_inbin = length(data$X)
num_part = data$particles.detected[1]
for (i in 1:length(data$X)) {
  num_inbin[i] = round(data$retained....[i] / 100 * num_part)
}
#jittering the data
n = sum(num_inbin)
counter = 0
jittered_data = rep(0, n)
data$endpoint[52] = 3.25
for (i in 1:length(num_inbin)) {
  interval_length = data$endpoint[i] - data$startpoint[i]
  for(j in 1:num_inbin[i]) {
    jittered_data[j + counter] = data$startpoint[i] + 
      runif(1, min = 0, max = interval_length)
  }
  counter = counter + num_inbin[i]
}
breaks = seq(from = 0, to = 3, by = 0.01)
hist(jittered_data, breaks=1000)
#testdata = rbilnorm(10000, 1, 5, 1, 1, 0.6)
#define starting parameters of the bilognormal model
mu1 <- -3
mu2 <- -0.5
sigma1 <- 0.3
sigma2 <- 0.4
tau <- 0.7
#estim = EM(testdata, mu1, mu2, sigma1, sigma2, tau)
#print(estim)
estim = EM(jittered_data, mu1, mu2, sigma1, sigma2, tau)
print(estim)


#optimization 
#write log likelihood function of binned data
cdfbilognorm = function(x, mu1, mu2, sigma1, sigma2, tau) {
  y <- (1 - tau) * plnorm(x, mu1, sigma1) + tau * plnorm(x, mu2, sigma2)
  return(y)
}

negloglikelihood = function(par) {
  logli = 1
  mu1 = par[1]
  mu2 = par[2]
  sigma1 = par[3]
  sigma2 = par[4]
  tau = par[5]
  for (i in 1:length(data$X)) {
    logli = logli * (cdfbilognorm(data$endpoint[i], mu1, mu2, sigma1, sigma2, tau) - 
                       cdfbilognorm(data$startpoint[i], mu1, mu2, sigma1, sigma2, tau)) ^ (num_inbin[i])
  }
  return(-log(logli))
}
out = optim(estim, negloglikelihood)


