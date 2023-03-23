rm(list = ls())


path_prefix = 'D:/HBI.Modeling'
project_name = 'MIMIC_III_ICU/modeling'
disease_name = 'sepsis'
version = 'version1'

source('D:/HBI.Modeling/ICU/helper_functions/func_global_visualization_helper_functions.r')
source('D:/HBI.Modeling/ICU/helper_functions/visualization.r')
source('D:/HBI.Modeling/ICU/helper_functions/general_helper.r')

rdata_path = paste(path_prefix, project_name, disease_name, version, 'rdata', sep = '/')
result_path = paste(path_prefix, project_name, disease_name, version, 'results', sep = '/')
delivery_path = paste(result_path, 'validation_results', sep = '/')
dir.create(delivery_path, showWarnings = F, recursive = T)

load(paste0(rdata_path, '/training_results.rdata'))
load(paste0(rdata_path, '/validation_data.rdata'))

validation_predicted_prob = predict(model, newdata = as.matrix(validation_matrix))
names(validation_predicted_prob) = rownames(validation_matrix)

#-----1. plot roc-----
global_visualization_helper_ROC_helper(validation_label, validation_predicted_prob, delivery_path, 'validation_ROC')

global_visualization_helper_boxplot(validation_predicted_prob, validation_label, label_name = c('control', 'case'), label_color = c('seagreen', 'violetred'),
                                    png_title = 'Validation Boxplot', delivery_path, png_name = 'validation_boxplot', if_remove_outliers = F)

save(validation_predicted_prob, file = paste0(rdata_path, '/validation_results.rdata'))




