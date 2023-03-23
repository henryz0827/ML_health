rm(list = ls())
library(xgboost)
library(ggplot2)


path_prefix = 'D:/HBI.Modeling/ICU/ICU_mimiciii'
disease_name = 'aki'
project_name = 'modeling'
version = 'version2'

rdata_path = paste(path_prefix, disease_name, project_name, version, 'rdata', sep = '/')
results_path =  paste(path_prefix, disease_name, project_name, version, 'results', sep = '/')
delivery_path = paste(results_path, 'training_results', sep = '/')

source('D:/HBI.Modeling/ICU/helper_functions/func_global_visualization_helper_functions.r')
source('D:/HBI.Modeling/ICU/helper_functions/visualization.r')
source('D:/HBI.Modeling/ICU/helper_functions/general_helper.r')


load(paste0(rdata_path, '/training_data.rdata'))
load(paste0(rdata_path, '/training_results.rdata'))
load(paste0(rdata_path, '/feature_info_total.rdata'))
group_weight_mat = data.frame(group = rownames(importance), importance = as.numeric(importance$Gain), stringsAsFactors = F)

group_weight_mat$group <- factor(group_weight_mat$group, levels = group_weight_mat$group[order(-group_weight_mat$importance)])

png(paste(delivery_path, '/feature_importance.png',sep = ''), width = 10, height = 8, unit = 'in', res = 400)
g = ggplot(data=group_weight_mat, aes(x = group, y = importance)) +
  geom_bar(stat="identity", fill="black")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(disease_name) +
  theme(plot.title = element_text(hjust = 0.5))
plot(g)
dev.off()

feature_info_total_append = feature_info_total
feature_info_total_append$used = 0
feature_info_total_append$used[rownames(feature_info_total_append) %in% group_weight_mat$group] = 1
write.csv(feature_info_total_append, file = paste0(results_path, '/feature_info_total_append.csv'))
save(feature_info_total_append, file = paste0(rdata_path, '/feature_info_total_append.rdata'))



