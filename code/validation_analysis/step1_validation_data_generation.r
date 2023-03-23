rm(list = ls())
library(xgboost)
library(ggplot2)


path_prefix = 'D:/HBI.Modeling/thuang/ICU/ICU_mimiciii'
disease_name = 'aki'
project_name = 'modeling'
version = 'version1'

rdata_path = paste(path_prefix, disease_name, project_name, version, 'rdata', sep = '/')
results_path =  paste(path_prefix, disease_name, project_name, version, 'results', sep = '/')
delivery_path = paste(results_path, 'validation_results', sep = '/')

load(paste0(rdata_path, '/data.rdata'))
load(paste0(rdata_path, '/validation_icu_samples.rdata'))
load(paste0(rdata_path, '/training_data.rdata'))
load(paste0(rdata_path, '/training_results.rdata'))

validation_sam_table = icu_sam_table[label$ICUSTAY_ID %in% validation_icu_samples, 
                                     match(colnames(training_matrix), colnames(icu_sam_table))]
validation_label = label[label$ICUSTAY_ID %in% validation_icu_samples, ]

save(validation_sam_table, validation_label, file = paste0(rdata_path, '/validation_data.rdata'))