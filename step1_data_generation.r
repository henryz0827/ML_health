rm(list = ls())
library(RODBC)

path_prefix = 'D:/HBI.Modeling/ICU/ICU_mimiciii'
disease_name = 'aki'
project_name = 'modeling'
version = 'version2'

rdata_path = paste(path_prefix, disease_name, project_name, version, 'rdata', sep = '/')
results_path =  paste(path_prefix, disease_name, project_name, version, 'results', sep = '/')
dir.create(rdata_path, recursive = T, showWarnings = F)
dir.create(results_path, recursive = T, showWarnings = F)

dbhandler <- odbcDriverConnect("driver={SQL Server};server=hbi-cndev-sql-1;
                               database=HBI_DATA_MIMIC_NEW;Trusted_Connection=yes")

#------0. save label icu_info------
load(paste(path_prefix, disease_name, 'generate_label/rdata/label.rdata', sep = '/'))
label$label = label$aki_label
save(label, file = paste0(rdata_path, "/label.rdata"))

#------0. remove icu sample with Y =1 in current time interval------
myfunc_case <- function(vec){
  if((length(vec) > 12) & (sum(vec[c(1:12)]) == 0) & (max(vec) == 1)){
    return(1)
  }
  return(0)
}

myfunc_control <- function(vec){
  if((length(vec) > 12) & (sum(vec[c(1:12)]) == 0) & (max(vec) == 0)){
    return(1)
  }
  return(0)
}

tmp_case = aggregate(label$label, by=list(label$ICUSTAY_ID),FUN= function(z) myfunc_case(z))
tmp_control = aggregate(label$label, by=list(label$ICUSTAY_ID),FUN= function(z) myfunc_control(z))
passed_sample_case = tmp_case$Group.1[tmp_case$x == 1]
passed_sample_control = tmp_control$Group.1[tmp_control$x == 1]
set.seed(12138)
used_sample = c(passed_sample_case, sample(passed_sample_control, length(passed_sample_case)))
save(used_sample, file = paste0(rdata_path, '/used_sample.rdata'))

label = label[label$ICUSTAY_ID %in% used_sample, ]
save(label, file = paste0(rdata_path, "/label_current_control_only.rdata"))

#------1. get icu_info------
used_icu = unique(label$ICUSTAY_ID)
cmd_icu_info  = paste0("SELECT [SUBJECT_ID],[ICUSTAY_ID],[INTIME],[OUTTIME],[LOS]
                        FROM [HBI_DATA_MIMIC_NEW].[Staging].[ICUSTAYS]
                       where [ICUSTAY_ID] in (", paste0(used_icu, collapse = ','), ")")
icu_info <- sqlQuery(dbhandler, cmd_icu_info, stringsAsFactors = F)
icu_info$INTIME = strptime(icu_info$INTIME, format = "%Y-%m-%d %H:%M:%S")
icu_info$OUTTIME = strptime(icu_info$OUTTIME, format = "%Y-%m-%d %H:%M:%S")
save(icu_info, file = paste0(rdata_path, '/icu_info.rdata'))

#------2. generate feature candidates ------ 
cmd_get_chart_items =  "SELECT [ITEMID],[LABEL],[ABBREVIATION],[CATEGORY]
                        FROM [HBI_DATA_MIMIC_NEW].[Staging].[D_ITEMS]
                        where ([CATEGORY] like '%lab%' or [CATEGORY] like '%Vital%' or [CATEGORY] like '%Respiratory%') and  [PARAM_TYPE] like '%numeric%'"
feature_table <- sqlQuery(dbhandler, cmd_get_chart_items, stringsAsFactors = F)
feature_table$feature_label = paste0("chart_", feature_table$ITEMID)
save(feature_table, file = paste0(rdata_path, '/feature_table.rdata'))

#------3. generate icu chart events table ------
cmd_get_charevents = paste0("SELECT [SUBJECT_ID]
                ,[ICUSTAY_ID]
                ,[ITEMID]
                ,[CHARTTIME]
                ,[VALUE]
                ,[VALUENUM]
                ,[VALUEUOM]
            FROM [HBI_DATA_MIMIC_NEW].[mimiciii].[chartevents_subitems]
            where [ICUSTAY_ID] in (", paste(used_icu, collapse = ','), ")")
icu_chartevents <- sqlQuery(dbhandler, cmd_get_charevents, stringsAsFactors = F)
icu_chartevents$CHARTTIME = strptime(icu_chartevents$CHARTTIME, format = "%Y-%m-%d %H:%M:%S")
save(icu_chartevents, file = paste0(rdata_path, '/icu_chartevents.rdata'))

#------4. generate sample data------
para_log_path = paste0(results_path, '/para_log2')
dir.create(para_log_path, showWarnings = F, recursive = T)
source('D:/HBI.Modeling/thuang/ICU/ICU_mimiciii/aki/helper_functions/helper_functions.r')
# icu_chartevents_with_label_list = chartevents_join_label(icu_info, icu_chartevents, label, feature_table, para_log_path)
# save(icu_chartevents_with_label_list, file = paste0(rdata_path, '/icu_chartevents_with_label_list.rdata'))
icu_sam_table = generate_sam_table(icu_chartevents, label, feature_table, para_log_path)
save(icu_sam_table, file = paste0(rdata_path, '/icu_sam_table.rdata'))


#------5. get age and gender------
dbhandler <- odbcDriverConnect("driver={SQL Server};server=hbi-cndev-sql-1;
                                 database=HBI_DATA_MIMIC_NEW;Trusted_Connection=yes")
cmd_subid_icuid =  paste0("SELECT [SUBJECT_ID], [ICUSTAY_ID] FROM [HBI_DATA_MIMIC_NEW].[Staging].[ICUSTAYS]
                        where [ICUSTAY_ID] in (", paste(used_icu, collapse = ","), ")")
subid_icuid_mapping <- sqlQuery(dbhandler, cmd_subid_icuid, stringsAsFactors = F)
all_subid = unique(subid_icuid_mapping$SUBJECT_ID)

cmd_subid_age_gender =  paste0("SELECT [SUBJECT_ID],[GENDER],[DOB]
                               FROM [HBI_DATA_MIMIC_NEW].[Staging].[PATIENTS]
                               where [SUBJECT_ID] in (", paste(all_subid, collapse = ","), ")")
subid_age_gender_mapping <- sqlQuery(dbhandler, cmd_subid_age_gender, stringsAsFactors = F)

age_gender_mapping = subid_age_gender_mapping[match(subid_icuid_mapping$SUBJECT_ID, subid_age_gender_mapping$SUBJECT_ID), c(2:3)]
subid_icuid_mapping = cbind(subid_icuid_mapping, age_gender_mapping)

label_mat_sufix = subid_icuid_mapping[match(label$ICUSTAY_ID, subid_icuid_mapping$ICUSTAY_ID), c(3, 4)]
label_mat_total = cbind(label,label_mat_sufix)

age = as.numeric(difftime(strptime(label_mat_total$Period_Start, format = "%Y-%m-%d %H:%M:%S"), 
         strptime(label_mat_total$DOB, format = "%Y-%m-%d %H:%M:%S"), unit = "weeks")) / 52.25
# male is 0, female is 1
gender = rep(0, nrow(label_mat_total))
gender[label_mat_total$GENDER == "F"] = 1

age_gender_df = data.frame(age = age, gender = gender)
icu_sam_table = cbind(age_gender_df, icu_sam_table)


save(icu_sam_table, label, 
     file = paste0(rdata_path, "/data.rdata"))



