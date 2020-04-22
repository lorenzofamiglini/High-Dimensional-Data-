library("caret")
library("tsutils")
library(ggplot2)

df <- read.csv("/Users/lorenzo/dev/uni/High-Dimensional-Data-/all_vars_df.csv")
str(df)
colnames(df)
head(df$id)
bmi_column <- 11
sex_column <-4
age_column <- 9
nutrients_columns <- 12:119
microbes_columns <- 120:238
pca_nutrients_columns <- 239:243
pcs_taxonomies_columns <- 244:249
cluster_column <- 250


# t test for sex regressor in bmi
x <-df[df$sex1m2f==1, 11]
y <-df[df$sex1m2f==2, 11]
t.test(x,y)


# names(getModelInfo())
#train.control <- trainControl(method = "cv", number = 3)
train.control <- trainControl(method = "cv", 
                              number = 10, savePredictions = "final", returnResamp = 'all')
print_lasso <- function(dataset) {
  # split in train and test: 90% obs in train and 10% in test
  train_n = floor(0.9 * nrow(dataset))
  test_n = floor(0.1 * nrow(dataset))
  train_ind = sample(seq_len(nrow(dataset)), size = train_n)
  training_set = dataset[train_ind,]
  test_set = dataset[-train_ind,]
  lambdas = lambdaseq(training_set[,-1], training_set[,1], weight = NA, alpha = 1, standardise = TRUE,
                            lambdaRatio = 1e-04, nLambda = 100, addZeroLambda = FALSE)$lambda 
  
  # Train the model
  model_lasso <- train(bmi ~., data = training_set, method = "glmnet", preProcess="scale",
                       trControl = train.control, tuneGrid=expand.grid(alpha=1, lambda=lambdas))
  #print(model_lasso)
  print(coef(model_lasso$finalModel,model_lasso$finalModel$lambdaOpt))
  mape = 100*mean(abs((test_set[,1] - as.numeric(predict(model_lasso, test_set)))/test_set[,1]))
  mse = mean((test_set[,1] - as.numeric(predict(model_lasso, test_set)))**2)
  
  resempla_cv <- model_lasso$resempla_cv <- model_lasso$resample
  mean_agg <- aggregate(resempla_cv$RMSE, list(resempla_cv$lambda), mean)
  std_cv <- aggregate(resempla_cv$RMSE, list(resempla_cv$lambda), sd)
  
  std_cv <- create_df(std_cv, met = "std")
  mean_cv <- create_df(mean_agg, met = "mean")
  mn_sd_cv <- cbind(mean_cv, std_cv)
  mn_sd_cv <- mn_sd_cv[,-3]
  
  mn_sd_cv$std_error <- 1/(sqrt(10))*mn_sd_cv$std
  mn_sd_cv$lower_bound <- mn_sd_cv$mean - 1.96 *mn_sd_cv$std_error
  mn_sd_cv$upper_bound <- mn_sd_cv$mean + 1.96 *mn_sd_cv$std_error
  
  res_r2 <- model_lasso$results$Rsquared
  
  mn_sd_cv$rsquared <- res_r2
  
  ci_stepwise_plot <- ggplot(mn_sd_cv) + 
    geom_point(data = mn_sd_cv, aes(x =`Regressor_#`,y = mean,color = "#000099"), alpha = 0.5) + 
    geom_errorbar(data = mn_sd_cv,aes(x = `Regressor_#`,y = mean, ymin = lower_bound, ymax = upper_bound,color = "#000099"), alpha = 0.5, width=.3,position=position_dodge(0.1))+
    geom_line(aes(x = `Regressor_#`,y = rsquared*10,color = "green")) + 
    geom_point(aes(x = `Regressor_#`,y = rsquared*10,color = "green")) + 
    scale_y_continuous(sec.axis = sec_axis(trans = ~./10,name = "R-squared")) + 
    labs(title="Stepwise Backward regression (CV 10 folds)",
         subtitle="Intervallo di confidenza CV rispetto al RMSE, R-squared (valor medio)",
         y = "RMSE (mean cv)",
         x = "Numero di regressori",
         color = "Metrica")+
    scale_color_manual(values=c("blue", "black"), labels = c("CI 95% RMSE", "R-squared"))
  
  
  
  
  lambdas_rmse = demo_nutrients_cluster_model$resample$RMSE
  mean_lambdas_rmse = mean(model_lasso$results$RMSE)
  sd_lambdas_rmse = sd(model_lasso$results$RMSE)
  list = list()
  list[["mape"]] = mape
  list[["mse"]] = mse
  list[["rmse"]] = mse**(1/2)
  list[["rmse_lower"]] = mn_sd_cv$lower_bound
  list[["rmse_upper"]] = mn_sd_cv$upper_bound
  print(list)
  return(ci_stepwise_plot)
}

create_df <- function(x, met = "std"){
  colnames(x)[1] <- "Regressor_#"
  colnames(x)[2] <- met
  return(x)
}

# print_ridge <- function(dataset) {
#   lambdas = lambdaseq(dataset[,-1], dataset[,1], weight = NA, alpha = 0, standardise = TRUE,
#                             lambdaRatio = 1e-04, nLambda = 100, addZeroLambda = FALSE)$lambda 
#   # Train the model
#   model_ridge <- train(bmi ~., data = dataset, method = "glmnet", preProcess="scale",
#                        trControl = train.control, tuneGrid=expand.grid(alpha=0, lambda=lambdas))
#   print(model_ridge)
#   print(coef(model_ridge$finalModel,model_ridge$finalModel$lambdaOpt ))
# }
# 
# 
# print_lm <- function(dataset) {
#   model_lm <- train(bmi ~., data = dataset, method = "lm", preProcess="scale",
#                     trControl = train.control)  
#   summary(model_lm)
# }


dataset <- df[, c(bmi_column, nutrients_columns, microbes_columns, sex_column, age_column, cluster_column)]
# togliamo zbmius e zbmicatus poiche' contenenti troppi na
summary(dataset[, c(2,3)])
trash_column = c(2,3,4)
dataset2 <- dataset[,-trash_column]
na_cord = which(is.na(dataset2), arr.ind = T)
dataset2[na_cord] = 0
# one hot encoding
clusters = names(table(dataset2[, 228]))
dataset2[dataset2$cluster==clusters[1], 229] = 0
dataset2[dataset2$cluster==clusters[2], 229] = 1
dataset3 = dataset2[,-228]
demo_taxa_nutr_cluster_model = print_lasso(dataset3)
write.csv(dataset3, file="demo_taxa_nutr_cluster.csv")
bmi_column2 = 1
nutrients_columns2 = 2:106
microbes_columns2 = 107:225
sex_column2 = 226
age_column2 = 227
cluster_column2 = 228
demo_taxa_cluster = dataset3[, c(bmi_column2, microbes_columns2, sex_column2, age_column2, cluster_column2 )]
demo_nutrients_cluster = dataset3[, c(bmi_column2, nutrients_columns2, sex_column2, age_column2, cluster_column2 )]
write.csv(demo_taxa_cluster, file="demo_taxa_cluster.csv")
write.csv(demo_nutrients_cluster, file="demo_nutrients_cluster.csv")
demo_taxa_cluster_model = print_lasso(demo_taxa_cluster)
demo_nutrients_cluster_model = print_lasso(demo_nutrients_cluster)


