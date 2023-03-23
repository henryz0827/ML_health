generate_label <- function(cur_icu_id, in_time, out_time, cur_time, icu_creatinine_info){
  # print(cur_time)
  if((cur_time + 48 * 3600) > out_time){
    return(NA)
  }
  
  start_time = cur_time
  end_time = cur_time + 48 * 3600
  
  sub_creatinine_info = icu_creatinine_info[icu_creatinine_info$ICUSTAY_ID == cur_icu_id, ]
  
  sub_creatinine_info = sub_creatinine_info[sub_creatinine_info$CHARTTIME >= start_time &
                                              sub_creatinine_info$CHARTTIME <= end_time, ]
  if(nrow(sub_creatinine_info) == 0){
    return(0)
  }
  sub_creatinine_info = sub_creatinine_info[!is.na(sub_creatinine_info$VALUE), 1:ncol(sub_creatinine_info), drop = F]
  if(nrow(sub_creatinine_info) == 0){
    return(0)
  }
  
  sub_creatinine_info = sub_creatinine_info[order(sub_creatinine_info$CHARTTIME), 1:ncol(sub_creatinine_info), drop = F]
  max_creatinine = max(sub_creatinine_info$VALUE)
  
  if(max_creatinine >= (1.5 * 1.2)){
    return(1)
  }
  
  if(nrow(sub_creatinine_info) == 1){
    return(0)
  }
  
  minimum = sub_creatinine_info$VALUE[1]
  diff = -1
  for(i in c(2:nrow(sub_creatinine_info))){
    cur_diff = sub_creatinine_info$VALUE[i] - minimum
    diff = max(diff, sub_creatinine_info$VALUE[i] - minimum)
    minimum = min(minimum, sub_creatinine_info$VALUE[i]);
  }
  if(diff >= 0.3){
    return(1)
  }
  return(0)
}



generate_label_forward <- function(cur_icu_id, in_time, out_time, cur_time, icu_creatinine_info){
  # print(cur_time)
  # if(cur_time <= in_time){
  #   return(NA)
  # }
  start_time = cur_time - 48 * 3600
  if(start_time < in_time){
    start_time = in_time
  }
  
  end_time = cur_time
  
  sub_creatinine_info = icu_creatinine_info[icu_creatinine_info$ICUSTAY_ID == cur_icu_id, ]
  
  sub_creatinine_info = sub_creatinine_info[sub_creatinine_info$CHARTTIME >= start_time &
                                              sub_creatinine_info$CHARTTIME <= end_time, ]
  if(nrow(sub_creatinine_info) == 0){
    return(0)
  }
  sub_creatinine_info = sub_creatinine_info[!is.na(sub_creatinine_info$VALUE), 1:ncol(sub_creatinine_info), drop = F]
  if(nrow(sub_creatinine_info) == 0){
    return(0)
  }
  
  sub_creatinine_info = sub_creatinine_info[order(sub_creatinine_info$CHARTTIME), 1:ncol(sub_creatinine_info), drop = F]
  max_creatinine = max(sub_creatinine_info$VALUE)
  
  if(max_creatinine >= (1.5 * 1.2)){
    return(1)
  }
  
  if(nrow(sub_creatinine_info) == 1){
    return(0)
  }
  
  minimum = sub_creatinine_info$VALUE[1]
  diff = -1
  for(i in c(2:nrow(sub_creatinine_info))){
    cur_diff = sub_creatinine_info$VALUE[i] - minimum
    diff = max(diff, sub_creatinine_info$VALUE[i] - minimum)
    minimum = min(minimum, sub_creatinine_info$VALUE[i]);
  }
  if(diff >= 0.3){
    return(1)
  }
  return(0)
}



generate_sam_table <- function(icu_chartevents, label, feature_table, para_log_path){
  
  icd_id = unique(label$ICUSTAY_ID)
  
  library(parallel)
  cl <- makeCluster(50)
  print("importing data")
  clusterExport(cl, c("icu_chartevents", "label", "feature_table", "icd_id", "para_log_path"), envir=environment())
  print("finish importing data")
  clusterEvalQ(cl, sink(paste0(para_log_path, '/', Sys.getpid(), ".txt")))
  icu_sam_table_list = parLapply(cl, 1:nrow(label), function(i){
    print(paste(i, '/', nrow(label)))
    icu_id = label$ICUSTAY_ID[i]
    Period_Start = label$Period_Start[i]
    Period_End = label$Period_End[i]
    
    sample_avg_vec = rep(NA, nrow(feature_table))
    sample_min_vec = rep(NA, nrow(feature_table))
    sample_max_vec = rep(NA, nrow(feature_table))
    
    cur_chartevents = icu_chartevents[(icu_chartevents$ICUSTAY_ID == icu_id), c(1:ncol(icu_chartevents)), drop = F]
  
    tmp_chartevents = cur_chartevents[(cur_chartevents$CHARTTIME >= Period_Start &
                                         cur_chartevents$CHARTTIME <= Period_End), c(1:ncol(cur_chartevents)), drop = F]
    tmp_chartevents = tmp_chartevents[!is.na(tmp_chartevents$VALUE), c("SUBJECT_ID", "ICUSTAY_ID", "ITEMID", "VALUE"), drop = F]
    
    if(nrow(tmp_chartevents) != 0){
      avg_mat = aggregate(tmp_chartevents,
                          by = list(tmp_chartevents$ITEMID),
                          FUN = mean)
      max_mat = aggregate(tmp_chartevents,
                          by = list(tmp_chartevents$ITEMID),
                          FUN = max)
      min_mat = aggregate(tmp_chartevents,
                          by = list(tmp_chartevents$ITEMID),
                          FUN = min)
      
      sample_avg_vec[match(avg_mat$ITEMID, feature_table$ITEMID)] = avg_mat$VALUE
      sample_max_vec[match(max_mat$ITEMID, feature_table$ITEMID)] = max_mat$VALUE
      sample_min_vec[match(min_mat$ITEMID, feature_table$ITEMID)] = min_mat$VALUE
    }
   
    
    names(sample_avg_vec) = paste0(feature_table$feature_label, "_avg")
    names(sample_min_vec) = paste0(feature_table$feature_label, "_min")
    names(sample_max_vec) = paste0(feature_table$feature_label, "_max")
    
    sample_vec = t(as.matrix(c(sample_min_vec, sample_max_vec, sample_avg_vec)))
      
    return(sample_vec)
  })
 
  stopCluster(cl)
  print("end paralel, start rbind")
  icu_sam_table = do.call(rbind, icu_sam_table_list)
  print("end rbind")
  return(icu_sam_table)
}


NA_imputation <- function(icu_sam_table, label, para_log_path){
  icu_sample = unique(label$ICUSTAY_ID)
  icu_sam_table = as.matrix(icu_sam_table)
  library(parallel)
  cl <- makeCluster(50)
  print("importing data")
  clusterExport(cl, c("icu_sam_table", "label", "icu_sample", "para_log_path"), envir=environment())
  print("finish importing data")
  clusterEvalQ(cl, sink(paste0(para_log_path, '/', Sys.getpid(), ".txt")))
  icu_sam_table_list = parLapply(cl, 1:length(icu_sample), function(i){
    print(paste(i, '/', length(icu_sample)))
    cur_icu = icu_sample[i]
    cur_sam_table = icu_sam_table[label$ICUSTAY_ID %in% cur_icu, ]
    for(j in 1:ncol(cur_sam_table)){
      cur_vec = cur_sam_table[, j]
      cur_vec_new = cur_vec
      if(sum(is.na(cur_vec)) != 0 & sum(is.na(cur_vec)) != length(cur_vec)){
        ind_value = c(0, which(!is.na(cur_vec)), length(cur_vec)+1)
        for(k in 1:(length(ind_value) - 1)){
          if((ind_value[k+1] - ind_value[k]) != 1){
            if(k == 1){
              cur_vec_new[1:(ind_value[k+1] - 1)] = cur_vec[ind_value[k+1]]
            }else if(k == length(ind_value) - 1){
              cur_vec_new[(ind_value[k] + 1):length(cur_vec_new)] = cur_vec[ind_value[k]]
            }else{
              mid_ind = floor((ind_value[k] + ind_value[k+1]) / 2)
              cur_vec_new[(ind_value[k] + 1):mid_ind] = cur_vec[ind_value[k]]
              cur_vec_new[(mid_ind + 1):(ind_value[k+1])] = cur_vec[ind_value[k+1]]
            }
          }
        }
      }
      cur_sam_table[, j] = cur_vec_new
    }
    return(cur_sam_table)
  })
  stopCluster(cl)
  print("end paralel, start rbind")
  icu_sam_table = do.call(rbind, icu_sam_table_list)
  print("end rbind")
  return(icu_sam_table)

}


feature_analysis_classification <- function(my_table, y){
  tmp_list = apply(my_table, 2, function(x){
    control_vec = x[y == 0]
    case_vec = x[y == 1]
    control_vec = control_vec[!is.na(control_vec)]
    case_vec = case_vec[!is.na(case_vec)]
    fc = median(case_vec) / median(control_vec)
    pval = wilcox.test(case_vec, control_vec)$p.value
    
    if(!is.na(fc) & fc < 1){
      tmp_pred = 1 / c(control_vec, case_vec)
    }else{
      tmp_pred = c(control_vec, case_vec)
    }
    pred = prediction(tmp_pred, c(rep(0, length(control_vec)), rep(1, length(case_vec))))
    perf <- performance(pred, "auc" )
    auc = round(unlist(perf@y.values),4)
    
    tmp_df = data.frame(fc = fc,
                        pval = pval,
                        auc = auc)
  })
  feature_info = do.call(rbind, tmp_list)
  return(feature_info)
}

chartevents_join_label <- function(icu_info, icu_chartevents, label, feature_table, para_log_path){
  icu_chartevents$INTIME = icu_info$INTIME[match(icu_chartevents$ICUSTAY_ID, icu_info$ICUSTAY_ID)]
  icu_chartevents$LOS = as.numeric(difftime(icu_chartevents$CHARTTIME, icu_chartevents$INTIME, units = "days"))
  
  label$INTIME = icu_info$INTIME[match(label$ICUSTAY_ID, icu_info$ICUSTAY_ID)]
  label$LOS_start = as.numeric(difftime(label$Period_Start, label$INTIME, units = "days"))
  label$LOS_end = as.numeric(difftime(label$Period_End, label$INTIME, units = "days"))
  
  #remove label sample 
  # cur_label = label[label$label == 0, ]
  unique_icu_sample = unique(label$ICUSTAY_ID)
  library(parallel)
  cl <- makeCluster(50)
  print("importing data")
  clusterExport(cl, c("icu_chartevents", "label", "feature_table", "unique_icu_sample", "para_log_path"), envir=environment())
  print("finish importing data")
  clusterEvalQ(cl, sink(paste0(para_log_path, '/', Sys.getpid(), ".txt")))
  icu_chartevents_with_label_list = parLapply(cl, 1:length(unique_icu_sample), function(i){
    print(paste0(i, '/', length(unique_icu_sample)))
    library("fuzzyjoin")
    cur_icu_sample = unique_icu_sample[i]
    sub_icu_chartevents = icu_chartevents[icu_chartevents$ICUSTAY_ID %in% cur_icu_sample, c("ICUSTAY_ID", "ITEMID", "VALUE", "LOS")]
    sub_label = label[label$ICUSTAY_ID %in% cur_icu_sample, c("ICUSTAY_ID", "Period_Start", "Period_End", 
                                                              "label", "LOS_start" , "LOS_end")]
    tmp = fuzzy_left_join(
      sub_icu_chartevents, sub_label,
      by = c(
        "LOS" = "LOS_start",
        "LOS" = "LOS_end"
      ),
      match_fun = list(`>=`, `<=`))
    tmp = as.matrix(tmp)
    return(tmp)
    
  })
  stopCluster(cl)
  names(icu_chartevents_with_label_list) = unique_icu_sample
  return(icu_chartevents_with_label_list)
  
}