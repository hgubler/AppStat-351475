load("2_online_shopping.RData")
library(dplyr)
library(ggplot2)
library(pROC)
data = Data # rename to "data"
rm(Data)
data = data %>% rename(purchase = Revenue, administrative = Administrative, administrative_duration = Administrative_Duration, informational = Informational, informational_duration = Informational_Duration, product_related = ProductRelated, product_related_duration = ProductRelated_Duration, bounce_rates = BounceRates, exit_rates = ExitRates, page_values = PageValues, special_day = SpecialDay, month = Month, operating_systems = OperatingSystems, browser = Browser, region = Region, traffic_type = TrafficType, visitor_type = VisitorType, weekend = Weekend) 
str(data)
#9 regions - is that too much to consider it as factor?
#20 traffic_types - is that too much to consider as factor?
#or should i just change the chr to factor variables and leave everything else as it is?
#change operating_systems, browser to factors
data$operating_systems = as.factor(data$operating_systems)
data$browser = as.factor(data$browser)
gm_all = glm(purchase ~ ., data = data)
res = gm_all$residuals
ggplot(data = data, aes(y = res, x = value)) + facet_wrap( ~ name, scales = "free") + 
  geom_point() + geom_smooth()
plot(gm_all)


AUC_eval <- function(gmodel,Data){
  set.seed(517)
  Folds <- matrix(sample(1:dim(Data)[1]), ncol=5)
  AUC <- rep(0,5)
  for(k in 1:5){
    train <- Data[-Folds[,k],]
    test <- Data[Folds[,k],]
    my_gm <- glm(gmodel$formula, family="binomial", data=train)
    test_pred <- predict(my_gm, newdata = test, type="response")
    AUC[k] <- auc(test$purchase,test_pred)
  }
  return(mean(AUC))
}
AUC_eval(gm_all, data)
