library("caret")
library("tsutils")

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


# t test for sex regressor in bmi
x <-df[df$sex1m2f==1, 11]
y <-df[df$sex1m2f==2, 11]
t.test(x,y)


# names(getModelInfo())
#train.control <- trainControl(method = "cv", number = 3)
train.control <- trainControl(method = "repeatedcv", 
                              number = 10, repeats = 3)
print_lasso <- function(dataset) {
  lambdas = lambdaseq(dataset[,-1], dataset[,1], weight = NA, alpha = 1, standardise = TRUE,
                            lambdaRatio = 1e-04, nLambda = 100, addZeroLambda = FALSE)$lambda 
  # Train the model
  model_lasso <- train(bmi ~., data = dataset, method = "glmnet", preProcess="scale",
                       trControl = train.control, tuneGrid=expand.grid(alpha=1, lambda=lambdas))
  print(model_lasso)
  print(coef(model_lasso$finalModel,model_lasso$finalModel$lambdaOpt ))
}

print_ridge <- function(dataset) {
  lambdas = lambdaseq(dataset[,-1], dataset[,1], weight = NA, alpha = 0, standardise = TRUE,
                            lambdaRatio = 1e-04, nLambda = 100, addZeroLambda = FALSE)$lambda 
  # Train the model
  model_ridge <- train(bmi ~., data = dataset, method = "glmnet", preProcess="scale",
                       trControl = train.control, tuneGrid=expand.grid(alpha=0, lambda=lambdas))
  print(model_ridge)
  print(coef(model_ridge$finalModel,model_ridge$finalModel$lambdaOpt ))
}


print_lm <- function(dataset) {
  model_lm <- train(bmi ~., data = dataset, method = "lm", preProcess="scale",
                    trControl = train.control)  
  summary(model_lm)
}

dataset <- df[,c(bmi_column, pca_nutrients_columns)]

print_lasso(dataset)
print_ridge(dataset)
print_lm(dataset)


dataset2 <- df[, c(bmi_column, nutrients_columns, sex_column, age_column)]
# non numeric
dataset2 <- dataset2[,-c(2,3)]
# replace na with 0
cord = which(is.na(dataset2[,c(1:99)]), arr.ind = T)
dataset2[cord] = 0
# remove likely collinear variable: bmicat1norm2ow3ob
dataset2 <- dataset2[, -2]
str(dataset2)
print_lasso(dataset2)
print_ridge(dataset2)
# lm non puo funzionare perche' ci sono troppi regressori
print_lm(dataset2)
