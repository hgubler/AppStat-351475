library("forecast")
library("tidyverse")
# load data and prepare it for time series
data_temp = read.csv(file = "data/temperature.csv")
data_temp = pivot_longer(data_temp, cols = c("Jan", "Feb", "Mar", 
                                             "Apr", "May", "Jun",
                                             "Jul", "Aug", "Sep", 
                                             "Oct", "Nov", "Dec"), 
                         names_to = "month", values_to = "temp")
data_temp <- subset(data_temp, Year >= 1950)
temp = as.numeric(data_temp$temp) # numerical vector of raw values
ts_temp = ts(temp, start = 1950, frequency = 12) # time series with monthly frequency
# calculate predictive checking error (1 step)
n = length(temp)
train = 1:(n-floor(n/10)) # training data (2 / 3 of real data)
err = array(0, c(5, floor(n/10)-1+1))
for(j in 0:(floor(n/10)-1)) {
  ts = ts(temp[train + j], start = 1950, frequency = 12)
  # fit potential models from analysis as well as suggestion from auto.arima
  ts_temp_fit_1_0_train = Arima(ts, order = c(2, 1, 0), seasonal = c(2, 0, 1), include.drift = TRUE)
  ts_temp_fit_2_0_train = Arima(ts, order = c(2, 2, 0), seasonal = c(2, 0, 1), include.drift = TRUE)
  ts_temp_fit_1_1_train = Arima(ts, order = c(2, 1, 0), seasonal = c(2, 1, 1), include.drift = TRUE)
  ts_temp_fit_2_1_train = Arima(ts, order = c(2, 2, 0), seasonal = c(2, 1, 1), include.drift = TRUE)
  ts_temp_fit_autoselect_train = Arima(ts, order = c(2, 1, 3), seasonal = c(1, 0, 0), include.drift = TRUE)
  err[1,j+1] = (temp[1 + end(train)[1]+j] - 
                  forecast(ts_temp_fit_1_0_train, h=1)$mean)^2
  err[2,j+1] = (temp[1 + end(train)[1]+j] - 
                  forecast(ts_temp_fit_2_0_train, h=1)$mean)^2
  err[3,j+1] = (temp[1 + end(train)[1]+j] - 
                  forecast(ts_temp_fit_1_1_train, h=1)$mean)^2
  err[4,j+1] = (temp[1 + end(train)[1]+j] - 
                  forecast(ts_temp_fit_2_1_train, h=1)$mean)^2
  err[5,j+1] = (temp[1 + end(train)[1]+j] - 
                  forecast(ts_temp_fit_autoselect_train, h=1)$mean)^2
  print(j)
}
mean_pred = as.data.frame(rowMeans(err))
save(mean_pred, file = "mean_pred.RData")
