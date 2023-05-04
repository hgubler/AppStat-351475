library("forecast")
# load data and prepare it for time series
data_temp = read.csv(file = "data/temperature.csv")
data_temp = pivot_longer(data_temp, cols = c("Jan", "Feb", "Mar", 
                                             "Apr", "May", "Jun",
                                             "Jul", "Aug", "Sep", 
                                             "Oct", "Nov", "Dec"), 
                         names_to = "month", values_to = "temp")
temp = as.numeric(data_temp$temp) # numerical vector of raw values
ts_temp = ts(temp, start = 1880, frequency = 12) # time series with monthly frequency
# fit potential models
ts_temp_fit_1_0 = Arima(ts_temp, order = c(2, 1, 0), seasonal = c(2, 0, 1))
ts_temp_fit_2_0 = Arima(ts_temp, order = c(2, 2, 0), seasonal = c(2, 0, 1))
ts_temp_fit_1_1 = Arima(ts_temp, order = c(2, 1, 0), seasonal = c(2, 1, 1))
ts_temp_fit_2_1 = Arima(ts_temp, order = c(2, 2, 0), seasonal = c(2, 1, 1))
# calculate predictive checking error (1 step)
n = length(temp)
train = 1:(n-floor(n/3)) # training data (2 / 3 of real data)
err = array(0, c(4, floor(n/3)-1+1))
for(j in 0:(floor(n/3)-1)) {
  ts_temp_fit_1_0_train = Arima(ts(temp[train + j], start = 1880, frequency = 12), 
                                order = c(2, 1, 0), seasonal = c(2, 0, 1))
  ts_temp_fit_2_0_train = Arima(ts(temp[train + j], start = 1880, frequency = 12), 
                                order = c(2, 2, 0), seasonal = c(2, 0, 1))
  ts_temp_fit_1_1_train = Arima(ts(temp[train + j], start = 1880, frequency = 12), 
                                order = c(2, 1, 0), seasonal = c(2, 1, 1))
  ts_temp_fit_2_1_train = Arima(ts(temp[train + j], start = 1880, frequency = 12), 
                                order = c(2, 2, 0), seasonal = c(2, 1, 1))
  err[1,j+1] = (temp[1 + end(train)[1]+j] - 
                  predict(ts_temp_fit_1_0_train, n.ahead=1)$pred)^2
  err[2,j+1] = (temp[1 + end(train)[1]+j] - 
                  predict(ts_temp_fit_2_0_train, n.ahead=1)$pred)^2
  err[3,j+1] = (temp[1 + end(train)[1]+j] - 
                  predict(ts_temp_fit_1_1_train, n.ahead=1)$pred)^2
  err[4,j+1] = (temp[1 + end(train)[1]+j] - 
                  predict(ts_temp_fit_2_1_train, n.ahead=1)$pred)^2
  print(j)
}