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
delivery_path = paste(result_path, 'training_results', sep = '/')
dir.create(delivery_path, showWarnings = F, recursive = T)

load(paste0(rdata_path, '/training_data.rdata'))
load(paste0(rdata_path, '/training_results.rdata'))

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


load(paste0(rdata_path, '/sub_feature_info.rdata'))
selected_feature_order = rownames(importance)[order(importance[, 1], decreasing = T)]

model_feature_info = sub_feature_info[match(selected_feature_order, make.names(rownames(sub_feature_info))), ]

write.csv(model_feature_info, file = paste0(delivery_path, '/model_feature_info.csv'))





