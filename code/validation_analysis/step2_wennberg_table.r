rm(list = ls())

path_prefix = 'D:/HBI.Modeling/tliu'
project_name = 'MIMIC_III_ICU/modeling'
disease_name = 'sepsis'
version = 'version1'

rdata_path = paste(path_prefix, project_name, disease_name, version, 'rdata', sep = '/')
result_path = paste(path_prefix, project_name, disease_name, version, 'results', sep = '/')
delivery_path = paste(result_path, 'validation_results', sep = '/')
dir.create(delivery_path, showWarnings = F, recursive = T)

source('D:/HBI.Modeling/tliu/ICU/helper_functions/general_helper.r')

load(paste0(rdata_path, '/training_data.rdata'))
load(paste0(rdata_path, '/training_results.rdata'))
load(paste0(rdata_path, '/validation_data.rdata'))
load(paste0(rdata_path, '/validation_results.rdata'))

# prevalence = 0.11

#-----2. bootstrap------
# risk_ratio = prevalence
# prob_label_list = bootstrap_population_risk(validation_predicted_prob, validation_label, risk_ratio)
# validation_predicted_prob_bootstrap = prob_label_list[[1]]
# validation_label_bootstrap = prob_label_list[[2]]

#---- wenberg table with bin ----
bin_vec = c(seq(0, 90, 10), 100.1) / 100
# prob = validation_predicted_prob
# label = validation_label_bootstrap
prob = validation_predicted_prob
label = validation_label
for(i in seq(1, length(bin_vec)-1)){
  # print(i)
  cur_bin = c(bin_vec[i], bin_vec[i + 1])
  
  medium_total_n = sum(prob >= cur_bin[1] & prob < cur_bin[2])
  medium_case = sum(label[prob >= cur_bin[1] & prob < cur_bin[2]] == 1)
  medium_control = sum(label[prob >= cur_bin[1] & prob < cur_bin[2]] == 0)
  
  medium_ppv = medium_case / medium_total_n
  medium_npv = medium_control / medium_total_n
  # medium_npv = (sum(label == 0) - medium_control) / (length(label) - medium_total_n)
  library(epiR)
  dat <- as.table(matrix(c(medium_case, sum(label == 0)-medium_control, sum(label == 1) - medium_case, medium_control), nrow = 2, byrow = TRUE))
  colnames(dat) <- c("Dis+","Dis-")
  rownames(dat) <- c("Test+","Test-")
  rval <- epi.tests(dat, conf.level = 0.95)
  sensitivity = round(rval$elements$sensitivity$est, 3)
  sensitivity_CI = paste0(round(rval$elements$sensitivity[, 2], 3), '~', round(rval$elements$sensitivity[, 3], 3))
  specificity = round(rval$elements$specificity$est, 3)
  specificity_CI = paste0(round(rval$elements$specificity[, 2], 3), '~', round(rval$elements$specificity[, 3], 3))
  
  if(i == 1){
    cutoff_table = data.frame(
      Intervals = c(paste0('[', cur_bin[1], ', ', cur_bin[2], ']')),
      Total_n = medium_total_n,
      Case_n = medium_case,
      PPV = medium_ppv,
      NPV = medium_npv,
      Sensitivity = sensitivity,
      Sensitivity_CI = sensitivity_CI,
      Specificity = specificity,
      Specificity_CI = specificity_CI,
      stringsAsFactors = F)
    
  }else{
    tmp = data.frame(
      Intervals = c(paste0('[', cur_bin[1], ', ', cur_bin[2], ']')),
      Total_n = medium_total_n,
      Case_n = medium_case,
      PPV = medium_ppv,
      NPV = medium_npv,
      Sensitivity = sensitivity,
      Sensitivity_CI = sensitivity_CI,
      Specificity = specificity,
      Specificity_CI = specificity_CI,
      stringsAsFactors = F)
    
    cutoff_table = rbind(cutoff_table, tmp)
  }
  
}

write.csv(cutoff_table, file = paste0(delivery_path, '/validation_wennberg.csv'))

# ------ wenberg table combine ----
low_medium_cutoff = 0.5
medium_high_cutoff = 0.5

bin_vec = c(0, low_medium_cutoff, medium_high_cutoff, 1.1)
prob = validation_predicted_prob
label = validation_label
for(i in seq(1, length(bin_vec)-1)){
  # print(i)
  cur_bin = c(bin_vec[i], bin_vec[i + 1])
  
  medium_total_n = sum(prob >= cur_bin[1] & prob < cur_bin[2])
  medium_case = sum(label[prob >= cur_bin[1] & prob < cur_bin[2]] == 1)
  medium_control = sum(label[prob >= cur_bin[1] & prob < cur_bin[2]] == 0)
  
  medium_ppv = medium_case / medium_total_n
  medium_npv = medium_control / medium_total_n
  
  library(epiR)
  dat <- as.table(matrix(c(medium_case, sum(label == 0)-medium_control, sum(label == 1) - medium_case, medium_control), nrow = 2, byrow = TRUE))
  colnames(dat) <- c("Dis+","Dis-")
  rownames(dat) <- c("Test+","Test-")
  rval <- epi.tests(dat, conf.level = 0.95)
  sensitivity = round(rval$elements$sensitivity$est, 3)
  sensitivity_CI = paste0(round(rval$elements$sensitivity[, 2], 3), '~', round(rval$elements$sensitivity[, 3], 3))
  specificity = round(rval$elements$specificity$est, 3)
  specificity_CI = paste0(round(rval$elements$specificity[, 2], 3), '~', round(rval$elements$specificity[, 3], 3))
  
  
  if(i == 1){
    cutoff_table_combine = data.frame(
      Intervals = c(paste0('[', cur_bin[1], ', ', cur_bin[2], ']')),
      Total_n = medium_total_n,
      Case_n = medium_case,
      PPV = medium_ppv,
      NPV = medium_npv,
      Sensitivity = sensitivity,
      Sensitivity_CI = sensitivity_CI,
      Specificity = specificity,
      Specificity_CI = specificity_CI,
      stringsAsFactors = F)
    
  }else{
    tmp = data.frame(
      Intervals = c(paste0('[', cur_bin[1], ', ', cur_bin[2], ']')),
      Total_n = medium_total_n,
      Case_n = medium_case,
      PPV = medium_ppv,
      NPV = medium_npv,
      Sensitivity = sensitivity,
      Sensitivity_CI = sensitivity_CI,
      Specificity = specificity,
      Specificity_CI = specificity_CI,
      stringsAsFactors = F)
    
    cutoff_table_combine = rbind(cutoff_table_combine, tmp)
  }
  
}

write.csv(cutoff_table_combine, file = paste0(delivery_path, '/validation_wennberg_combine.csv'))


final_wennberg = generate_final_wennberg(validation_predicted_prob, validation_label, cutoff = low_medium_cutoff)
rownames(final_wennberg) = disease_name
write.csv(final_wennberg, file = paste0(delivery_path, '/final_wennberg.csv'))

# ----- normal distribution------ 
set.seed(12138)
validation_predicted_prob_normal_distribution = rnorm(length(validation_predicted_prob), 
                                                             mean(validation_predicted_prob), 
                                                             sd(validation_predicted_prob))
tmp_prevalence = sum(validation_predicted_prob >= medium_high_cutoff) / length(validation_predicted_prob)
tmp_prevalence2 = sum(validation_predicted_prob >= low_medium_cutoff) / length(validation_predicted_prob)
low_medium_cutoff_normal_distribution = quantile(validation_predicted_prob_normal_distribution, probs = c(1-tmp_prevalence2))
medium_high_cutoff_normal_distribution = quantile(validation_predicted_prob_normal_distribution, probs = c(1-tmp_prevalence))

# --- save data----
save(validation_predicted_prob_normal_distribution,
     file = paste0(rdata_path, '/validation_results_normal_distribution.rdata'))

save(low_medium_cutoff, medium_high_cutoff, low_medium_cutoff_normal_distribution, medium_high_cutoff_normal_distribution,
     file = paste0(rdata_path, '/validation_cutoff.rdata'))

