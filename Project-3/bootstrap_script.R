library(MASS)
load("Project-3/Data/bootstrap_data.RData")
gm_full_inter = glm(fulltime_goals ~ . + home:covid, data = data, family = "poisson")
set.seed(27) # for reproducability
n_boot = 1000 # number of bootstrap
fisher_mat = vcov(gm_full_inter) # get the inverse of the fisher information
model_mat = model.matrix(gm_full_inter) # model matrix
deviance_boot = rep(0, n_boot) # to store all the deviances of the bootstrapped models
# sample the bootstrapped coefeficients
coef_sample = mvrnorm(n = n_boot, mu = unname(coef(gm_full_inter)), Sigma = fisher_mat)
for (i in 1:n_boot) {
  # simulate response data using poisson distibution and means from the bootstrapped
  # models
  goals_boot = rpois(nrow(data), lambda = exp(model_mat %*% coef_sample[i,]))
  # create bootstrapped dataset with sampled goals from above
  data_boot = data
  data_boot$fulltime_goals = goals_boot
  # poisson regression model with poisson sampled goals as response from above
  gm_boot = glm(fulltime_goals ~ . + home:covid, data = data_boot, 
                family = "poisson")
  # store deviance for i-th bootstrapped model
  deviance_boot[i] = gm_boot$deviance
}
# calculate p-value estimate
p = length(deviance_boot[deviance_boot > (gm_full_inter$deviance)]) / n_boot

# export the p_value as a dataframe
# p_value = as.data.frame(p)
# save(p_value, file = "bootstrap_data.RData")

