rm(list = ls())

path_prefix = 'D:/HBI.Modeling/tliu'
project_name = 'MIMIC_III_ICU/modeling'
disease_name = 'sepsis'
version = 'version1'

source('D:/HBI.Modeling/tliu/ICU/helper_functions/func_global_visualization_helper_functions.r')
source('D:/HBI.Modeling/tliu/ICU/helper_functions/visualization.r')
source('D:/HBI.Modeling/tliu/ICU/helper_functions/general_helper.r')

rdata_path = paste(path_prefix, project_name, disease_name, version, 'rdata', sep = '/')
result_path = paste(path_prefix, project_name, disease_name, version, 'results', sep = '/')
delivery_path = paste(result_path, 'validation_results', sep = '/')


load(paste0(rdata_path, '/training_data.rdata'))
load(paste0(rdata_path, '/training_results.rdata'))
load(paste0(rdata_path, '/validation_data.rdata'))
load(paste0(rdata_path, '/validation_results.rdata'))


png_path = paste0(delivery_path, '/ppv_plot.png')
ppv_score_df = ppv_plot(png_path, title_name = 'score vs ppv', validation_predicted_prob, validation_label, 
                    fit_model = 'loess', risk_ratio = NA, if_smooth = F)
save(ppv_score_df, file = paste0(rdata_path, '/validation_ppv_score_df.rdata'))










