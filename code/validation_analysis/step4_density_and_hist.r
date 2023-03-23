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


load(paste0(rdata_path, '/training_results.rdata'))
load(paste0(rdata_path, '/validation_data.rdata'))
load(paste0(rdata_path, '/validation_results.rdata'))
load(paste0(rdata_path, '/validation_cutoff.rdata'))
load(paste0(rdata_path, '/validation_results_normal_distribution.rdata'))

global_visualization_helper_single_density_plot(validation_predicted_prob, testing_sample_val = NA,
                                                delivery_path, cur_trend = 'up', cur_cp = disease_name, if_has_testing = F, if_remove_outlier = F, 
                                                png_name = 'density', is_vertical_line = T, 
                                                vertical_line_threshold1 = low_medium_cutoff, vertical_line_threshold2 = medium_high_cutoff,
                                                high_color = 'violetred', low_color = 'seagreen', medium_color = 'yellow')

global_visualization_helper_single_density_plot(validation_predicted_prob_normal_distribution, testing_sample_val = NA,
                                                delivery_path, cur_trend = 'up', cur_cp = disease_name, if_has_testing = F, if_remove_outlier = F, 
                                                png_name = 'density_normal_distribution', is_vertical_line = T, 
                                                vertical_line_threshold1 = low_medium_cutoff_normal_distribution, vertical_line_threshold2 = medium_high_cutoff_normal_distribution,
                                                high_color = 'violetred', low_color = 'seagreen', medium_color = 'yellow')

global_visualization_helper_single_hist_plot(testing_sample_val = NA, prevalence = NA,
                                             delivery_path, cur_trend = 'up', cur_cp = disease_name, if_has_testing = F, png_name = 'hist', 
                                             low_medium_cutoff, medium_high_cutoff, 
                                             validation_predicted_prob, baseline_label = NA, high_color = 'violetred', low_color = 'seagreen', medium_color = 'yellow', 
                                             cutoff_change = T, if_no_medium = T)

global_visualization_helper_two_hist_with_density_plot(validation_predicted_prob[validation_label == 0], 
                                                       validation_predicted_prob[validation_label == 1],
                                                       delivery_path, bin_color = 'darkgoldenrod1', 
                                                       xlab = 'Predicted Probability', ylab = 'Frequency', 
                                                       title = paste0('Histogram of ', disease_name, ','),
                                                       num_breaks = 10, title_size = 2, lab_size = 2)

global_visualization_helper_scatter_plot(validation_predicted_prob, validation_label, 
                                        delivery_path, low_medium_cutoff)
