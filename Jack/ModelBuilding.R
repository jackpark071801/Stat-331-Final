setwd("D:/Waterloo/Waterloo 2B/Stat 331/Final Project")
pollutants <- read.csv("pollutants.csv")
pollutants$X <- NULL
factors <- c(27, 28, 29, 32)
for (i in 1:length(factors)) {
  pollutants[,factors[i]] <- as.factor(pollutants[,factors[i]])
}

library(car)
model <- lm(length ~ ., pollutants)
vif(model)
maxVIF <- 10
tempdata <- pollutants
while (maxVIF >= 10) {
  model <- lm(length ~ ., tempdata)
  maxVIF <- max(vif(model))
  if (maxVIF >= 10) {
    tempdata <- tempdata[,- (which.max(vif(model))+1)]
  }
}
data_VIF <- tempdata

model <- lm(length ~ ., data_VIF)
summary(model)

set.seed(331)
index <- sample(1:nrow(data_VIF), 0.7*nrow(data_VIF))
train <- data_VIF[index,]
test <- data_VIF[-index,]
M <- lm(length ~ ., train)
M_test <- lm(length ~ ., test)

eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
}

library(glmnet)
lambdas <- 10^seq(2, -3, by = -0.1)
x <- model.matrix(M)
y_train <- train$length
lasso_reg <- cv.glmnet(x, y_train, alpha = 1, lambda = lambdas)
lambda_best <- lasso_reg$lambda.min

lasso_model <- glmnet(x, y_train, alpha = 1, lambda = lambda_best, standardize = TRUE)
predictions_train <- predict(lasso_model, s = lambda_best, newx = x)
eval_results(y_train, predictions_train, train)

x_test <- model.matrix(M_test)
y_test <- test$length
predictions_test <- predict(lasso_model, s = lambda_best, newx = x_test)
eval_results(y_test, predictions_test, test)

lasso_model$beta

plot(y_test, predictions_test, main = "LASSO Model", xlab = "Length",
     ylab = "Predicted Length")

#k fold cross validation, p-value, AIC and BIC

library(caret)
train.control <- trainControl(method = "LOOCV")
kFoldModel <- train(length ~., data = pollutants, method = "lm",
                    trControl = train.control)
print(kFoldModel)
