library("knitr")
library("car")
library("nlme")
# load data and prepare it for time series
data_temp = read.csv(file = "data/temperature.csv")
data_temp = pivot_longer(data_temp, cols = c("Jan", "Feb", "Mar", 
                                             "Apr", "May", "Jun",
                                             "Jul", "Aug", "Sep", 
                                             "Oct", "Nov", "Dec"), 
                         names_to = "month", values_to = "temp")
# only take observations later than 1950
data_temp <- subset(data_temp, Year >= 1950)
# create dataframe for gls
lm_data = data.frame(y = temp, 
                     month = as.factor(rep(1:12,length(temp)/12)),
                     time = 1:length(temp))
# calculate correlattion structure using p and q parameters from sarima model
corr_struct = corARMA(form = ~ time, p = 2, q = 3)
# fit gls and calculate anova for square dependence
glsfit_lin = gls(y ~ time + month, data = lm_data, corr = corr_struct, method = "ML")
lin_anova = Anova(glsfit_lin, type = 2)
# fit gls and calculate anova for square dependence
glsfit_square = gls(y ~ time + I(time^2) + month, data = lm_data, corr = corr_struct, method = "ML")
square_anova = Anova(glsfit_square, type = 2)
# fit gls and calclulate anova for exp dependence
lm_data$time = lm_data$time / 10 # so we don't have too big values, otherwise get inf for exp(time)
glsfit_exp = gls(y ~ time + I(exp(time)) + month, data = lm_data, corr = corr_struct, method = "ML")
exp_anova = Anova(glsfit_exp, type = 2)
lin_anova = as.data.frame(lin_anova)
square_anova = as.data.frame(square_anova)
exp_anova = as.data.frame(exp_anova)
save(lin_anova, file = "lin_anova.RData")
save(square_anova, file = "square_anova.RData")
save(exp_anova, file = "exp_anova.RData")
