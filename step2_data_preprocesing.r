rm(list = ls())
library(ROCR)

path_prefix = 'D:/HBI.Modeling/ICU/ICU_mimiciii'
disease_name = 'aki'
project_name = 'modeling'
version = 'version2'

rdata_path = paste(path_prefix, disease_name, project_name, version, 'rdata', sep = '/')
results_path =  paste(path_prefix, disease_name, project_name, version, 'results', sep = '/')
dir.create(rdata_path, recursive = T, showWarnings = F)
dir.create(results_path, recursive = T, showWarnings = F)

source('D:/HBI.Modeling/ICU/ICU_mimiciii/aki/helper_functions/helper_functions.r')

time_interval = 48 * 2
load(paste0(rdata_path, '/data.rdata'))

#--------0. matrix imputation--------
para_log_path = paste0(results_path, '/log_imputation')
dir.create(para_log_path, showWarnings = F, recursive = T)
icu_sam_table_imputed = NA_imputation(icu_sam_table, label, para_log_path)
icu_sam_table_imputed = data.frame(icu_sam_table_imputed)
save(icu_sam_table_imputed, file = paste0(rdata_path, '/icu_sam_table_imputed.rdata'))

#--------1. get qualified sample(current is control)--------
actual_sam_table = icu_sam_table[c(1: (nrow(label) - time_interval)), ]
actual_label = label[c(1: (nrow(label) - time_interval)), ]
shift_label = label[c((time_interval + 1):nrow(label)), ]
used_sample_ind = which((actual_label$ICUSTAY_ID == shift_label$ICUSTAY_ID) &
                          (difftime(shift_label$Period_Start, actual_label$Period_Start,units = "days") == (time_interval / 24)))

sub_sam_table = actual_sam_table[used_sample_ind, ]
sub_shift_label = shift_label[used_sample_ind, ]
sub_actual_label = actual_label[used_sample_ind, ]

sub_used_row = which(sub_actual_label$label == 0)
sub_sam_table = sub_sam_table[sub_used_row, ]
sub_shift_label = sub_shift_label[sub_used_row, ]
sub_actual_label = sub_actual_label[sub_used_row, ]

#--------2. sample feature na remove--------
final_sam_table = sub_sam_table[, colSums(is.na(sub_sam_table)) / nrow(sub_sam_table) < 0.5]
tmp_ind = which(rowSums(is.na(final_sam_table)) / ncol(final_sam_table) < 0.5)
final_sam_table = final_sam_table[tmp_ind, ]
final_actual_label = sub_actual_label[tmp_ind, ]
final_shift_label = sub_shift_label[tmp_ind, ]

#--------3. feature analysis--------
feature_info = feature_analysis_classification(final_sam_table, final_shift_label$label)
save(feature_info, file = paste0(rdata_path, '/feature_info.rdata'))

#--------4. separate train / validation icu sample--------
case_icu_sample = unique(final_shift_label$ICUSTAY_ID[final_shift_label$label == 1])
control_icu_sample = setdiff(unique(final_shift_label$ICUSTAY_ID), case_icu_sample)

set.seed(12138)
train_case_sample = sample(case_icu_sample, length(case_icu_sample) / 2)
set.seed(12138)
train_control_sample = sample(control_icu_sample, length(control_icu_sample) / 2)

train_ind = which(final_shift_label$ICUSTAY_ID %in% c(train_case_sample, train_control_sample))
training_matrix = final_sam_table[train_ind, ]
training_actual_label =  final_actual_label$label[train_ind]
training_shift_label = final_shift_label$label[train_ind]
training_icu_info = cbind(final_actual_label[train_ind, ], final_shift_label[train_ind, ])

validation_matrix = final_sam_table[-train_ind, ]
validation_actual_label =  final_actual_label$label[-train_ind]
validation_shift_label = final_shift_label$label[-train_ind]
validation_icu_info = cbind(final_actual_label[-train_ind, ], final_shift_label[-train_ind, ])

library(impute)
training_matrix_imputed = impute.knn(as.matrix(training_matrix),k = 10)$data
feature_info_training = feature_analysis_classification(training_matrix_imputed, training_shift_label)
save(feature_info_training, file = paste0(rdata_path, '/feature_info_training.rdata'))

save(training_matrix_imputed, training_matrix, training_actual_label, training_shift_label, training_icu_info,
     file = paste0(rdata_path, '/training_data.rdata'))
save(validation_matrix, validation_actual_label, validation_shift_label,validation_icu_info,
     file = paste0(rdata_path, '/validation_data.rdata'))

#--------5. append feature info---------
load(paste0(rdata_path, '/feature_table.rdata'))
feature_more_info = feature_table[match(substr(rownames(feature_info), 1, nchar(rownames(feature_info))-4), feature_table$feature_label), ]
feature_info_total = cbind(feature_info, feature_more_info)
feature_info_total = feature_info_total[order(feature_info_total$auc, decreasing = T), ]

write.csv(feature_info_total, file = paste0(results_path, '/feature_info_total.csv'))
save(feature_info_total, file = paste0(rdata_path, '/feature_info_total.rdata'))
