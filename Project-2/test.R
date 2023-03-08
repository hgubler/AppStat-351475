load("Project-2/2_online_shopping.RData")
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(forcats)
library(splines)
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

data = Data # rename to "data"
rm(Data)
data = data %>% rename(purchase = Revenue, administrative = Administrative,
                       administrative_duration = Administrative_Duration, 
                       informational = Informational, informational_duration = Informational_Duration,
                       product_related = ProductRelated, product_related_duration = ProductRelated_Duration,
                       bounce_rates = BounceRates, exit_rates = ExitRates,
                       page_values = PageValues, special_day = SpecialDay, month = Month,
                       operating_systems = OperatingSystems, browser = Browser, region = Region,
                       traffic_type = TrafficType, visitor_type = VisitorType, weekend = Weekend) 
str(data)
data = data %>% mutate_if(is.character, as.factor)
data = data %>% mutate_if(is.integer, as.numeric)
data = data %>% mutate_if(is.logical, as.factor)
data = data %>% mutate_at(c("visitor_type", "region", "browser", "operating_systems",
                            "traffic_type"), as.factor) 

#levels(data$month) = c(levels(data$month), "other")
#tt = table(data$browser)
#data$browser = fct_collapse(data$browser, "other" = names(tt[tt < 40]))
for (i in 1:length(data)) {
  if(is.factor(data[, i])) {
    tt = table(data[, i])
    data[, i] = fct_collapse(data[, i], "other" = names(tt[tt < 40]))
  }
}

gm_all = glm(purchase ~ ., data=data, family="binomial")
# do facetwrap
data %>% mutate(res = resid(gm_all)) %>% mutate_if(is.factor, as.numeric) %>% 
  mutate_if(is.integer, as.numeric) %>% pivot_longer(-res) %>% 
  ggplot(aes(y = res, x = value)) + facet_wrap( ~ name, scales = "free") + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ ns(x, 2))

AUC_eval(gm_all, data)

