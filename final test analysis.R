pollutants <- read.csv("C:/Users/Shawn/Downloads/pollutants.csv")
library("car")
library(ggplot2)
library(ggcorrplot)
library(caret)
library(glmnet)
#remove the x
pollutants["X"] = NULL

#calculate correlation matrix before removing the multicolinearity covariates
corr_mat = cor(pollutants)
#graph colored corr matrix
ggcorrplot(corr_mat)

#set the factor
pollutants$edu_cat = as.factor(pollutants$edu_cat)
pollutants$race_cat = as.factor(pollutants$race_cat)
pollutants$male = as.factor(pollutants$male)
pollutants$smokenow = as.factor(pollutants$smokenow)
#fit model
model = lm(length ~ .,data = pollutants)
#summary
summary(model)
#show the VIF
vif(model)

#get set a dataset with no categorical covariates
no_cat = pollutants
no_cat$edu_cat = NULL
no_cat$race_cat = NULL
no_cat$male = NULL
no_cat$smokenow = NULL
#summary of the dataset
summary(no_cat)
#calculate correlation matrix
corr_matrix = cor(no_cat)
#graph colored corr matrix
ggcorrplot(corr_matrix)

#show the qqPlotogram of all covariates
qqPlot(pollutants$length)
qqPlot(pollutants$POP_PCB1)
qqPlot(pollutants$POP_PCB2)
qqPlot(pollutants$POP_PCB3)
qqPlot(pollutants$POP_PCB4)
qqPlot(pollutants$POP_PCB5)
qqPlot(pollutants$POP_PCB6)
qqPlot(pollutants$POP_PCB7)
qqPlot(pollutants$POP_PCB8)
qqPlot(pollutants$POP_PCB9)
qqPlot(pollutants$POP_PCB10)
qqPlot(pollutants$POP_PCB11)
qqPlot(pollutants$POP_dioxin1)
qqPlot(pollutants$POP_dioxin2)
qqPlot(pollutants$POP_dioxin3)
qqPlot(pollutants$POP_furan1)
qqPlot(pollutants$POP_furan2)
qqPlot(pollutants$POP_furan3)
qqPlot(pollutants$POP_furan4)
qqPlot(pollutants$whitecell_count)
qqPlot(pollutants$lymphocyte_pct)
qqPlot(pollutants$monocyte_pct)
qqPlot(pollutants$eosinophils_pct)
qqPlot(pollutants$basophils_pct)
qqPlot(pollutants$neutrophils_pct)
qqPlot(pollutants$BMI)
qqPlot(pollutants$ageyrs)
qqPlot(pollutants$yrssmoke)
qqPlot(pollutants$ln_lbxcot)

#find the linearity (residual y against residual x)
covariates = names(no_cat)
for (name in covariates){
  y_model = lm(paste("length", "~", "."), data = pollutants)
  x_model = lm(paste(name, "~", ".", "- length"), data = pollutants)
  y_resid = resid(y_model)
  x_resid = resid(x_model)
  #linearity
  scatterplot(x_resid, y_resid, ylab = "Length", xlab = name)
}

#find the model fit for homoscedasticity before remove multicolinearity
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))

#remove covariates with VIF > 10 below:

#remove highest vif ~ 15047 
pollutants$eosinophils_pct = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#summary
summary(model)
#show the VIF
vif(model)

#remove highest vif eosinophils_pct ~ 15047 then POP_PCB5 ~ 59.718357
pollutants$POP_PCB5 = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#summary
summary(model)
#show the VIF
vif(model)

#remove highest vif eosinophils_pct ~ 15047 then POP_PCB5 ~ 59.718357 then POP_PCB1 ~ 32.255513
pollutants$POP_PCB1 = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#summary
summary(model)
#show the VIF
vif(model)

#remove highest vif eosinophils_pct ~ 15047 then POP_PCB5 ~ 59.718357 then POP_PCB1 ~ 32.255513 then POP_PCB2 ~ 22.564616
pollutants$POP_PCB2 = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#summary
summary(model)
#show the VIF
vif(model)

#remove highest vif eosinophils_pct ~ 15047 then POP_PCB5 ~ 59.718357 then POP_PCB1 ~ 32.255513 then POP_PCB2 ~ 22.564616
# then POP_PCB4 ~ 14.260201
pollutants$POP_PCB4 = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#summary
summary(model)
#show the VIF
vif(model)

#get set a dataset with no categorical covariates
no_cat = pollutants
no_cat$edu_cat = NULL
no_cat$race_cat = NULL
no_cat$male = NULL
no_cat$smokenow = NULL
#summary of the dataset
summary(no_cat)
#calculate correlation matrix
corr_matrix = cor(no_cat)
#graph colored corr matrix
ggcorrplot(corr_matrix)


#find the model fit for homoscedasticity after removing multicolinearity
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))

#set up train and test model
data = model.matrix(length ~ ., data = pollutants)
n = nrow(data)
#set seed for sample and 
set.seed(331)
#get index for random train and test index (90% train, 10% test)
train_row = sample(1:n, 0.9*n)
#train set
x_matrix = data[,-1]
y_matrix = pollutants$length

#get the y values
train_y = y_matrix[train_row]
#get the x values
train_x = x_matrix[train_row,]
#get the y values
test_y = y_matrix[-train_row]
#get the x values
test_x = x_matrix[-train_row,]

eval_results <- function(true, predicted, df){
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square)
}

#k fold lasso
lambdas <- 10^seq(2, -4, by = -.0001)
lasso = cv.glmnet(train_x, train_y, alpha = 1, lambda = lambdas)
best_lam = lasso$lambda.min
best_lam
#use lasso with best lambda
best_lasso = glmnet(train_x, train_y, alpha = 1, lambda = best_lam)
#predict with lasso result with training set
pred_train = predict(best_lasso, s = best_lam, newx = train_x)
eval_results(train_y, pred_train, train_x)

#predict with lasso result with test set
pred_train = predict(best_lasso, s = best_lam, newx = test_x)
eval_results(test_y, pred_train, test_x)

#training data set
train_data = pollutants[train_row,]
#training model
train_model = lm(length ~ ., data = train_data)

#step wise AIC on training data
step_aic = step(train_model, direction = "both")
summary(step_aic)
aic_pred = predict(step_aic, newdata = pollutants[-train_row,])

#RMSE AIC
aic_true = pollutants$length[-train_row]
aic_sd = sum((aic_true - aic_pred)^2)
msd_aic = aic_sd / length(aic_true)
rmse_aic = sqrt(msd_aic)
rmse_aic

#step wise BIC on training data
step_bic = step(train_model, direction = "both", k = log(nrow(train_x)))
summary(step_bic)
bic_pred = predict(step_bic, newdata = pollutants[-train_row,])

#RMSE AIC
bic_true = pollutants$length[-train_row]
bic_sd = sum((bic_true - bic_pred)^2)
msd_bic = bic_sd / length(bic_true)
rmse_bic = sqrt(msd_bic)
rmse_bic