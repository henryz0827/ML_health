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


load(paste0(rdata_path, '/training_data.rdata'))
load(paste0(rdata_path, '/training_results.rdata'))
load(paste0(rdata_path, '/validation_data.rdata'))
time_interval = 96
icu_samples = unique(validation_label$ICUSTAY_ID)
validation_sam_table_imputed = validation_sam_table

for(j in 1:ncol(validation_sam_table_imputed)){
  imputed_Value = median(training_matrix[training_label == 0, 
                                         match(colnames(validation_sam_table_imputed)[j], colnames(training_matrix))])
  validation_sam_table_imputed[is.na(validation_sam_table_imputed[, j]), j] = imputed_Value
}
validation_results = predict(model, as.matrix(validation_sam_table_imputed))


validation_label$predict = validation_results


#------population level future 48 hr, boxplot-----
prob_df = validation_label[1:(nrow(validation_label) - 1*time_interval), ]
label_df = validation_label[(1+1*time_interval):nrow(validation_label), ]
shared_ind = which(prob_df$ICUSTAY_ID == label_df$ICUSTAY_ID)
df = data.frame(prob_vec = as.numeric(prob_df$predict[shared_ind]),
                label = as.numeric(label_df$aki_label[shared_ind]))
df$label = as.factor(df$label)


png_folder = paste(delivery_path, 'population', sep = '/')
dir.create(png_folder, recursive = T, showWarnings = F)
png(paste0(png_folder, "/time_interval_", time_interval, ".png"))
p<-ggplot(df, aes(x=label, y=prob_vec, fill=label)) +
  geom_boxplot()
plot(p)
dev.off()


#------population level future 48 hr, roc-----
library(ROCR)
pred = prediction(df$prob_vec, df$label)
perf <- performance(pred, "auc" )
auc = round(unlist(perf@y.values), 4)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
png_folder = paste(delivery_path, 'population', sep = '/')
dir.create(png_folder, recursive = T, showWarnings = F)
png(paste0(png_folder, "/roc_time_interval_", time_interval, ".png"),
    width = 10, height = 10, units = 'in', res = 400)
plot(roc.perf, main = paste0('ROC ', round(auc, 2)), cex.main = 2, cex.lab = 1.3)
dev.off()

perff=performance(pred, "tpr", "fpr" )
min.idx <- which(((1-perff@x.values[[1]])^2 + (perff@y.values[[1]])^2)==max((1-perff@x.values[[1]])^2 + (perff@y.values[[1]])^2))
cutoff=pred@cutoffs[[1]][min.idx[1]]


#------population level future 48 hr, cutoff find-----
low_medium_cutoff = 0.4
medium_high_cutoff = 0.4

bin_vec = c(0, low_medium_cutoff, medium_high_cutoff, 1)
prob = df$prob_vec
label = df$label
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




#------sample level future 48 hr line-----
png_folder = paste(delivery_path, 'continuous_prob', sep = '/')
dir.create(png_folder, recursive = T, showWarnings = F)
for(i in 1:length(icu_samples)){
  print(paste(i, '/', length(icu_samples)))
  cur_icu = icu_samples[i]
  result_df = validation_label[validation_label$ICUSTAY_ID %in% cur_icu, ]
  result_df$Idx = as.numeric(result_df$Idx)

  label_df = result_df[(1+1*time_interval):nrow(result_df), ]
  
  
  if((1+1*time_interval) < nrow(result_df)){
    result_df$future_label = NA
    result_df$future_label[match(label_df$Idx - time_interval, result_df$Idx)] = label_df$aki_label
    png(paste0(png_folder, "/icuid_", i, "_", cur_icu, ".png"))
    p = ggplot(data=result_df)+ ylim(0, 1)+
      geom_line(mapping=aes(y=aki_label,x= Idx,color="label"),size=1 ) +
      geom_line(mapping=aes(y=predict,x= Idx,color="probability"),size=1) +
      geom_line(mapping=aes(y=future_label,x= Idx,color="future label"),size=1, linetype = 2) +
      geom_hline(yintercept=low_medium_cutoff, linetype="dashed", 
                 color = "brown", size=2)+
      scale_color_manual(values = c(
        'label' = 'darkblue',
        'probability' = 'red',
        'future label' = 'seagreen')) +
      xlab('icu stay time (hours)') + ylab("label and probability") +
      ggtitle(paste(disease_name, "\n icu id: ", cur_icu)) +
      theme(plot.title = element_text(hjust = 0.5)) 
    plot(p)
    dev.off()
  }else{
    png(paste0(png_folder, "/icuid_", i, "_", cur_icu, ".png"))
    p = ggplot(data=result_df)+ ylim(0, 1)+
      geom_line(mapping=aes(y=aki_label,x= Idx,color="label"),size=1 ) +
      geom_line(mapping=aes(y=predict,x= Idx,color="probability"),size=1) +
      geom_hline(yintercept=low_medium_cutoff, linetype="dashed", 
                 color = "brown", size=2)+
      scale_color_manual(values = c(
        'label' = 'darkblue',
        'probability' = 'red')) +
      xlab('icu stay time (hours)') + ylab("label and probability") +
      ggtitle(paste(disease_name, "\n icu id: ", cur_icu)) +
      theme(plot.title = element_text(hjust = 0.5)) 
    plot(p)
    dev.off()
  }
  
  
  
}

