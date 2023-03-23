rm(list = ls())
library(xgboost)
library(ggplot2)


path_prefix = 'D:/HBI.Modeling/thuang/ICU/ICU_mimiciii'
disease_name = 'aki'
project_name = 'modeling'
version = 'version2'

source('D:/HBI.Modeling/thuang/ICU/helper_functions/func_global_visualization_helper_functions.r')


rdata_path = paste(path_prefix, disease_name, project_name, version, 'rdata', sep = '/')
results_path =  paste(path_prefix, disease_name, project_name, version, 'results', sep = '/')
delivery_path = paste(results_path, 'training_results', sep = '/')
dir.create(delivery_path, showWarnings = F, recursive = T)

load(paste0(rdata_path, '/training_data.rdata'))
training_label = training_shift_label
d = xgb.DMatrix(data = as.matrix(training_matrix), label = as.numeric(as.character(training_label)))

max_error = 0
parameters_vec = c()
for(max_depth in 1:10){
  for(min_child_weight in 1:10){
    for(gamma in seq(0, 1, 0.1)){
      set.seed(12138)
      cv <- xgb.cv(data = d, nrounds = 100, nfold = 5, eval_metric = "auc", nthread = 5,
                   max_depth = max_depth,  
                   min_child_weight = min_child_weight,
                   gamma = gamma,
                   objective = "binary:logistic", verbose = 0, early_stopping_rounds = 50)
      best_error = as.numeric(cv$evaluation_log[cv$best_iteration, 4])
      print(paste(max_depth, min_child_weight, gamma,best_error, max_error))
      if(best_error > max_error){
        max_error = best_error
        print(max_error)
        parameters_vec = c(max_depth, min_child_weight, gamma)
      }
    }
  }
}

parameters_vec
max_depth = 6
min_child_weight = 1
gamma = 0
scale_pos_weight = round(sum(training_label == 0) / sum(training_label == 1), 2)
max_error = 0
for(eta in seq(0.1, 0.8, 0.1)){
  for(subsample in seq(0.5, 1, 0.1)){
    for(colsample_bytree in seq(0.5, 1, 0.1)){
      # for(scale_pos_weight in c(seq(0.1, 0.9, 0.1), seq(1, 10, 1))){
      set.seed(12138)
      cv <- xgb.cv(data = d, nrounds = 100, nfold = 5, eval_metric = "auc", nthread = 5,
                   eta = eta,
                   max_depth = max_depth,
                   min_child_weight = min_child_weight,
                   gamma = gamma,
                   scale_pos_weight = scale_pos_weight,
                   subsample = subsample,
                   colsample_bytree = colsample_bytree,
                   # booster = "gblinear",
                   objective = "binary:logistic", verbose = 0, early_stopping_rounds = 50)
      best_error = as.numeric(cv$evaluation_log[cv$best_iteration, 4])
      best_iteration = cv$best_iteration
      print(paste(eta, subsample, colsample_bytree, scale_pos_weight, best_error, max_error, best_iteration, parameters_vec[8]))
      if(best_error > max_error){
        max_error = best_error
        print(max_error)
        parameters_vec = c(eta, max_depth, min_child_weight, gamma,
                           subsample, colsample_bytree, scale_pos_weight, best_iteration)
      }
      
      # }
      
    }
    
  }
}

parameters_vec = c(  0.10,  2.00,  1.00 , 0.00,  0.50 , 0.70,  5.53, 60)

set.seed(12138)
cv <- xgb.cv(data = d, nrounds = parameters_vec[8], nfold = 5, eval_metric = "auc", nthread = 5,
             eta = parameters_vec[1],
             max_depth = parameters_vec[2],  
             min_child_weight = parameters_vec[3],
             gamma = parameters_vec[4],
             subsample = parameters_vec[5],
             colsample_bytree = parameters_vec[6],
             scale_pos_weight = parameters_vec[7],
             objective = "binary:logistic", verbose = 1, prediction = T)

set.seed(12138)
model = xgb.train(data = d, nrounds = parameters_vec[8], nthread = 5,
                  eta = parameters_vec[1],
                  max_depth = parameters_vec[2],  
                  min_child_weight = parameters_vec[3],
                  gamma = parameters_vec[4],
                  subsample = parameters_vec[5],
                  colsample_bytree = parameters_vec[6],
                  scale_pos_weight = parameters_vec[7],
                  objective = "binary:logistic")
training_predicted_prob = predict(model, as.matrix(training_matrix))

global_visualization_helper_boxplot(training_predicted_prob, training_label, label_name = c('control', 'case'), label_color = c('green', 'red'),
                                    png_title = 'training_boxplot', delivery_path, png_name = 'training_boxplot', if_remove_outliers = F)

importance = data.frame(xgb.importance(model = model))
rownames_importance = importance$Feature
importance = importance[, -1]
rownames(importance) = rownames_importance

global_visualization_helper_ROC_helper(training_label, training_predicted_prob, delivery_path, 'training_ROC')

save(model, importance, training_predicted_prob, file = paste0(rdata_path, '/training_results.rdata'))




