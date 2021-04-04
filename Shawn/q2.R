energyapp <- read.csv("C:/Users/Shawn/Downloads/energyapp.csv")
pollutants <- read.csv("C:/Users/Shawn/Downloads/pollutants.csv")

backward_adj_r2 = function(dataset, response, exclude = "")
{
  # Specify the columns that are exclude/eliminated from model after ~
  exclude = c(exclude, response)
  #get all the name of columns in the dataset
  cols = names(dataset)
  #create avaliable covariate list (covariates tht are not elimiated)
  col_list = c()
  #return model at the end
  true_model = NULL
  #filling the col_list (avaliable covariate list)
  for (col in cols)
  {
    #if the covariate is not eliminated/excluded
    if (!(col %in% exclude))
    {
      #add this covariate to the avaliable covariate list
      col_list = c(col_list, col)
    }
  }
  #all the covariates + reponse
  all_cols = c(cols,"")
  #set the max adjusted R^2 for all the possible model in this dataset
  max_adj_r2 = 0
  #set the previous adjusted R^2 value
  old_adj_r2 = -1
  #place holder for model
  model = lm(paste(response, "~", "."), data = dataset)
  #if there are still larger adjusted R^2 value when a covariate is removed
  while (max_adj_r2 > old_adj_r2)
  {
    #set the name placeholder for covariate being tested for removal
    selected_col = ""
    #keep the previous largest adjusted R^2
    old_adj_r2 = max_adj_r2
    #calculate adjusted R^2 value for each covariate in the model when removed
    for (col in all_cols)
    {
      #don't test on covariate already eliminated
      if (!(col %in% exclude))
      {
        #find all the covariates that are not eliminatied and is not col
        sel_cols = setdiff(col_list, col)
        #set all above covariates to linear regression equation
        sel_cols = paste(sel_cols,collapse = "+")
        #fit the model with response ~ sel_cols with dataset
        model = lm(formula = paste(response, " ~ ",sel_cols) , data = dataset)
        #get the summary of the model and find the adjusted R^2
        adj_r2 = summary(model)$adj.r.squared
        #if this model has the biggest adjusted R^2
        if (adj_r2 > max_adj_r2)
        {
          #set this to be the max adjusted R^2 so far
          max_adj_r2 = adj_r2
          #get the name of coefficient that gives the max adjusted R^2 value when removed form model
          selected_col = col
          true_model = model
        }
        #print(s$coefficients)
        col = if(col=="") "NONE" else col
        #write the covariate removed from model and the adjusted R^2 value after fitting the removed model.
        cat(format(paste("-",col), width = 15), " Adj-R^2:", format(adj_r2,width = 15, justify = 'right'),"\n")
      }
    }
    #set the name of covariate eliminated
    selected_col = if(selected_col=="") "NONE" else selected_col
    #print the elimination message along with the max adjusted R^2 value
    cat("\n==> Removed ",selected_col,": Adj-R^2:",max_adj_r2, "\n\n")
    #update the available covariate list
    col_list = setdiff(col_list, selected_col)
    #adding the eliminated covariate to the elimination list
    exclude = c(exclude, selected_col)
  }
  #return the final model
  return(true_model)
}


backward_p_value = function(dataset, response, alpha) {
  #creating a dataset to manipulate covariates
  varset = dataset
  #removing response in the covariates dataset
  varset[response] = NULL
  #creating the df for response use in lm
  response = dataset[response]
  
  #place holder for model
  model = NULL
  #loop while there are still p-values > alpha
  while (TRUE)
  {
    #max p-value for each lm call that is greater than alpha
    maxp = 0
    #creating lm for, use this model to delete a p-value > alpha
    model = lm(paste(response, "~", "."), data = varset)
    cat("\n\n", format(paste("New Model"), width = 15, justify = 'centre'), "\n")
    #getting the covariates of p-value
    pvalues = summary(model)$coefficients
    #index for the covariate to be deleted
    index = 0
    #number of covariates in the model
    ncols = ncol(varset)
    #checking all the p-values except for intercept
    for(ind in 2:(1+ncols)){
      #getting the p-values for each covariate
      pval = pvalues[ind, "Pr(>|t|)"]
      #checking if the p-value is greater than alpha and is a new max p-value
      if((pval > alpha) & (pval > maxp)) {
        #set this to be the max p-value
        maxp = pval
        #get the index, use to delete this covariate column later
        index = ind - 1
        cat(format(paste("-",names(varset[index])), width = 15), "p-value > alpha:", format(pval,width = 15, justify = 'right'),"\n")
      }
    }
    #when there are no p-values greater than alpha
    if(maxp <= 0){
      #return the model
      cat("\n\n", format(paste("Final Model Found!"), width = 15, justify = 'centre'), "\n")
      if (!is.null(model)){
        return(model)
      } else {
        return(lm(paste(response, "~", "."), data = dataset))
      }
    } else {
      #remove covariate that has the greatest p-value greater than alpha
      #deleting the column
      cat(format(paste("-",names(varset[index])), width = 15), " largest p-value:", format(pvalues[index+1, "Pr(>|t|)"],width = 15, justify = 'right'),"\n")
      varset[index] = NULL
    }
  }
}


backward_aic = function(dataset, response, exclude = ""){
  # Specify the columns that are exclude/eliminated from model after ~
  exclude = c(exclude, response)
  #get all the name of columns in the dataset
  cols = names(dataset)
  #create available covariate list (covariates tht are not eliminated)
  col_list = c()
  #filling the col_list (available covariate list)
  for (col in cols){
    #if the covariate is not eliminated/excluded
    if (!(col %in% exclude)){
      #add this covariate to the available covariate list
      col_list = c(col_list, col)
    }
  }
  #all the covariates + response
  all_cols = c(cols,"")
  #set the min AIC for all the possible model in this dataset
  full_model = lm(paste(response, "~", "."), data = dataset)
  min_aic = AIC(full_model)
  #set the previous AIC
  old_aic = min_aic + 1
  #place holder for model
  return_model = full_model
  #if there are still smaller AIC value when a covariate is removed
  while (min_aic < old_aic){
    #set the name placeholder for covariate being tested for removal
    selected_col = ""
    #keep the previous smallest AIC
    old_aic = min_aic
    comb = paste(col_list,collapse = "+")
    org_model = lm(paste(response, " ~ ",comb), data = dataset)
    cat(format(paste("=","Model"), width = 15), " AIC:", format(old_aic,width = 15, justify = 'right'),"\n")
    min_model = NULL
    #calculate AIC value for each covariate in the model when removed
    for (col in all_cols){
      #don't test on covariate already eliminated
      if (!(col %in% exclude)){
        #find all the covariates that are not eliminatied and is not col
        sel_cols = setdiff(col_list, col)
        #set all above covariates to linear regression equation
        sel_cols = paste(sel_cols,collapse = "+")
        #fit the model with response ~ sel_cols with dataset
        model = lm(formula = paste(response, " ~ ",sel_cols) , data = dataset)
        #get the summary of the model and find the AIC
        aic = AIC(model)
        #if this model has the smallest AIC
        if(aic < min_aic){
          #set this to be the min AIC so far
          min_aic = aic
          #get the name of coefficient that gives the min AIC value when removed form model
          selected_col = col
          min_model = model
        }
        col = if(col=="") "NONE" else col
        #write the covariate removed from model and the AIC value after fitting the removed model.
        if (aic < old_aic){
          cat(format(paste("-",col), width = 15), " AIC:", format(aic,width = 15, justify = 'right'),"\n")
        } else {
          cat(format(paste("+",col), width = 15), " AIC:", format(aic,width = 15, justify = 'right'),"\n")
        }
      }
    }
    #set the name of covariate eliminated
    selected_col = if(selected_col=="") "NONE" else selected_col
    #print the elimination message along with the min AIC value
    cat("\n==> Removed ",selected_col,": AIC:",min_aic, "\n\n")
    #update the avaliable covariate list
    col_list = setdiff(col_list, selected_col)
    #adding the eliminated covariate to the elimination list
    exclude = c(exclude, selected_col)
  }
  if(!is.null(min_model)) {
    return_model = min_model
  } else {
    return_model = org_model
  }
  #return the final model
  return(return_model)
}


backward_bic = function(dataset, response, exclude = ""){
  # Specify the columns that are exclude/eliminated from model after ~
  exclude = c(exclude, response)
  #get all the name of columns in the dataset
  cols = names(dataset)
  #create available covariate list (covariates tht are not eliminated)
  col_list = c()
  #filling the col_list (available covariate list)
  for (col in cols){
    #if the covariate is not eliminated/excluded
    if (!(col %in% exclude)){
      #add this covariate to the available covariate list
      col_list = c(col_list, col)
    }
  }
  #all the covariates + reponse
  all_cols = c(cols,"")
  #set the min BIC for all the possible model in this dataset
  full_model = lm(paste(response, "~", "."), data = dataset)
  min_bic = BIC(full_model)
  #set the previous BIC
  old_bic = min_bic + 1
  #place holder for model
  return_model = full_model
  
  #if there are still smaller BIC value when a covariate is removed
  while (min_bic < old_bic){
    #set the name placeholder for covariate being tested for removal
    selected_col = ""
    #keep the previous smallest BIC
    old_bic = min_bic
    cat(format(paste("=","Model"), width = 15), " BIC:", format(old_bic,width = 15, justify = 'right'),"\n")
    comb = paste(col_list,collapse = "+")
    org_model = lm(paste(response, " ~ ",comb), data = dataset)
    min_model = NULL
    #calculate BIC value for each covariate in the model when removed
    for (col in all_cols){
      #don't test on covariate already eliminated
      if (!(col %in% exclude)){
        #find all the covariates that are not eliminated and is not col
        sel_cols = setdiff(col_list, col)
        #set above covariates to linear regression equation
        sel_cols = paste(sel_cols,collapse = "+")
        #fit the model with response ~ sel_cols with dataset
        model = lm(formula = paste(response, " ~ ",sel_cols) , data = dataset)
        #get the summary of the model and find the BIC
        bic = BIC(model)
        #if this model has the smallest BIC
        if(bic < min_bic){
          #set this to be the min BIC so far
          min_bic = bic
          #get the name of coefficient that gives the min BIC value when removed form model
          selected_col = col
          min_model = model
        }
        #print(s$coefficients)
        col = if(col=="") "NONE" else col
        #write the covariate removed from model and the BIC value after fitting the removed model.
        if (bic < old_bic){
          cat(format(paste("-",col), width = 15), " BIC:", format(bic,width = 15, justify = 'right'),"\n")
        } else {
          cat(format(paste("+",col), width = 15), " BIC:", format(bic,width = 15, justify = 'right'),"\n")
        }
      }
    }
    if(!is.null(min_model)) {
      return_model = min_model
    } else {
      return_model = org_model
    }
    #set the name of covariate eliminated
    selected_col = if(selected_col=="") "NONE" else selected_col
    #print the elimination message along with the min BIC value
    cat("\n==> Removed ",selected_col,": BIC:",min_bic, "\n\n")
    #update the avaliable covariate list
    col_list = setdiff(col_list, selected_col)
    #adding the eliminated covariate to the elimination list
    exclude = c(exclude, selected_col)
  }
  #return the final model
  return(return_model)
}

# take dataset as the dataset
# response as string of the response variable
# method as "adjr2" as the adjusted R^2 method
# method as "AIC" as the AIC method
# method as "BIC" as the BIC method
# method as "pvalue" as the p-value method
# alpha as the alpha value for method p value, takes a number
# pvalue is not necessary for adjusted R^2, AIC or BIC method, only needed for p-value
# performs backward elimination on the dataset
# if NULL is return, then the full model is the best model according to backward elimination result
# Example: backward(energyapp, "appEuse", "adjr2")
#          backward(energyapp, "appEuse", "pvalue", 0.05)
backward = function(dataset, response, method, alpha) {
  #method using adjusted R^2
  if (method == "adjr2") {
    backward_adj_r2(dataset = dataset, response = response)
    #method using AIC
  } else if (method == "AIC") {
    backward_aic(dataset = dataset, response = response)
    #method using BIC
  } else if (method == "BIC") {
    backward_bic(dataset = dataset, response = response)
    #method suing pvalue
  } else if (method == "pvalue") {
    backward_p_value(dataset = dataset, response = response, alpha = alpha)
  }
}

#dataset from reduced model in Q3b
#create new model with VIF > 10 removed
energyapp[c("T_out", "T1", "T6", "RH_7", "RH_2", "T2", "Tdewpoint", "RH_4", "RH_9")] = NULL

#setting the index for the first 500 training data
index = seq(1,500)
#getting the first 500 training data
train = energyapp[index,]
#getting all the data after the first 500 for testing
test = energyapp[-index,]
test_row = 479

pollutants$X = NULL
#turn categorical numbers into strings
pollutants$male[pollutants$male == 1] = "Male"
pollutants$male[pollutants$male == 0] = "Female"

pollutants$smokenow[pollutants$smokenow == 1] = "Smokes"
pollutants$smokenow[pollutants$smokenow == 0] = "No Smokes"

pollutants$edu_cat[pollutants$edu_cat == 1] = "12 grade no dip below"
pollutants$edu_cat[pollutants$edu_cat == 2] = "High School Grad/GED,Eqv"
pollutants$edu_cat[pollutants$edu_cat == 3] = "Some College,AA degree"
pollutants$edu_cat[pollutants$edu_cat == 4] = "College graduate"

pollutants$race_cat[pollutants$race_cat == 1] = "Other race, mutli-racial"
pollutants$race_cat[pollutants$race_cat == 2] = "Mexican American"
pollutants$race_cat[pollutants$race_cat == 3] = "Non-Hispanic Black"
pollutants$race_cat[pollutants$race_cat == 4] = "Non-Hispanic White"

#factor the categorical
pollutants$male = as.factor(pollutants$male)
pollutants$smokenow = as.factor(pollutants$smokenow)
pollutants$edu_cat = as.factor(pollutants$edu_cat)
pollutants$race_cat = as.factor(pollutants$race_cat)

model = lm(length ~ .,data = pollutants)
#show the VIF
vif(model)


#remove covariates with VIF > 10 below:

#remove highest vif ~ 15047 
pollutants$eosinophils_pct = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#show the VIF
vif(model)

#remove highest vif eosinophils_pct ~ 15047 then POP_PCB5 ~ 59.718357
pollutants$POP_PCB5 = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#show the VIF
vif(model)

#remove highest vif eosinophils_pct ~ 15047 then POP_PCB5 ~ 59.718357 then POP_PCB1 ~ 32.255513
pollutants$POP_PCB1 = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#show the VIF
vif(model)

#remove highest vif eosinophils_pct ~ 15047 then POP_PCB5 ~ 59.718357 then POP_PCB1 ~ 32.255513 then POP_PCB2 ~ 22.564616
pollutants$POP_PCB2 = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#show the VIF
vif(model)

#remove highest vif eosinophils_pct ~ 15047 then POP_PCB5 ~ 59.718357 then POP_PCB1 ~ 32.255513 then POP_PCB2 ~ 22.564616
# then POP_PCB4 ~ 14.260201
pollutants$POP_PCB4 = NULL
#fit new model
model = lm(length ~ .,data = pollutants)
#show the VIF
vif(model)

#adjusted R^2
result1 = backward(pollutants, "length", method = "adjr2")
#summary for adjusted R^2 final model
summary(result1)


#AIC
result2 = backward(pollutants, "length", method = "AIC")
#summary for AIC final model
summary(result2)


#BIC
result3 = backward(pollutants, "length", method = "BIC")
#summary for BIC final model
summary(result3)


#p-value with alpha of 0.05
result4 = backward(pollutants, "length", method = "pvalue", alpha = 0.05)
#summary for p-value with alpha = 0.05 final model
summary(result4)
#AIC of the final model
