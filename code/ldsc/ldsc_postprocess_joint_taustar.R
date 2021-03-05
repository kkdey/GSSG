

###########################  many annots tau star  ##################################

library(data.table)
library(rmeta)
annot_cell = "/n/groups/price/kushal/Enhancer_MasterReg/data/ANNOTATIONS/FigJoints/Joint_All_genes"
results_cell = "/n/groups/price/kushal/Enhancer_MasterReg/data/LDSC_RESULTS/FigJoints/Joint_All_genes/baselineLD_v2.1_no_Coding_no_Promoter_no_TSS"
base_path = "/n/groups/price/kushal/Enhancer_MasterReg/data/ANNOTATIONS/Baselines/baselineLD_v2.1_no_Coding_no_Promoter_no_TSS"
annot_names = "FS1"


all_traits = c('UKB_460K.body_BMIz','UKB_460K.cov_EDU_YEARS','UKB_460K.lung_FVCzSMOKE','UKB_460K.cov_SMOKING_STATUS',
               'UKB_460K.mental_NEUROTICISM','UKB_460K.blood_WHITE_COUNT','PASS_Years_of_Education2','UKB_460K.bp_SYSTOLICadjMEDz',
               'UKB_460K.body_HEIGHTz','UKB_460K.other_MORNINGPERSON','UKB_460K.body_WHRadjBMIz','UKB_460K.lung_FEV1FVCzSMOKE',
               'UKB_460K.repro_MENARCHE_AGE','UKB_460K.blood_RED_COUNT','UKB_460K.blood_PLATELET_COUNT','UKB_460K.bmd_HEEL_TSCOREz',
               'UKB_460K.blood_EOSINOPHIL_COUNT','PASS_Schizophrenia','UKB_460K.blood_RBC_DISTRIB_WIDTH','PASS_Height1','PASS_BMI1',
               'UKB_460K.disease_T2D','PASS_AgeFirstBirth','UKB_460K.disease_RESPIRATORY_ENT','UKB_460K.body_BALDING1','UKB_460K.disease_HYPOTHYROIDISM_SELF_REP',
               'UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED','UKB_460K.disease_HI_CHOL_SELF_REP','UKB_460K.repro_MENOPAUSE_AGE','PASS_HDL','UKB_460K.pigment_SUNBURN',
               'PASS_NumberChildrenEverBorn','PASS_Anorexia','PASS_LDL','PASS_Crohns_Disease','PASS_DS','PASS_Ever_Smoked','UKB_460K.pigment_HAIR',
               'PASS_Rheumatoid_Arthritis','PASS_Type_2_Diabetes','PASS_Autism','UKB_460K.pigment_TANNING','PASS_Ulcerative_Colitis',
               'UKB_460K.disease_DERMATOLOGY','PASS_Coronary_Artery_Disease','UKB_460K.disease_AID_SURE','UKB_460K.pigment_SKIN')

blood_traits = c("UKB_460K.blood_RBC_DISTRIB_WIDTH", "UKB_460K.blood_RED_COUNT", "UKB_460K.blood_WHITE_COUNT",
                 "UKB_460K.blood_PLATELET_COUNT", "UKB_460K.blood_EOSINOPHIL_COUNT")

autoimmune_traits = c("UKB_460K.disease_AID_SURE", "PASS_Ulcerative_Colitis", "PASS_Crohns_Disease",
                      "PASS_Rheumatoid_Arthritis", "PASS_Celiac",  "PASS_Lupus", "PASS_Type_1_Diabetes",
                      "PASS_IBD", "PASS_Primary_biliary_cirrhosis")
# autoimmune_traits = c("PASS_IBD", "PASS_Ulcerative_Colitis", "PASS_Crohns_Disease",
#                       "PASS_Rheumatoid_Arthritis", "PASS_Celiac",  "PASS_Type_1_Diabetes")

brain_traits = c("PASS_Ever_Smoked", "UKB_460K.cov_SMOKING_STATUS", "UKB_460K.mental_NEUROTICISM", "UKB_460K.repro_MENARCHE_AGE",
                 "PASS_Years_of_Education2", "PASS_DS", "PASS_Schizophrenia", "UKB_460K.body_WHRadjBMIz",
                 "PASS_BMI1", "UKB_460K.body_BMIz")


immune_traits = c("PASS_Celiac", "PASS_Crohns_Disease", "PASS_IBD", "PASS_Lupus",
                  "PASS_Primary_biliary_cirrhosis", "PASS_Rheumatoid_Arthritis",
                  "PASS_Type_1_Diabetes", "PASS_Ulcerative_Colitis",
                  "UKB_460K.disease_ASTHMA_DIAGNOSED", "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
                  "PASS_Multiple_sclerosis", "UKB_460K.body_BMIz", "UKB_460K.disease_HYPERTENSION_DIAGNOSED",
                  "PASS_Triglycerides", "PASS_LDL", "PASS_HDL",
                  "UKB_460K.bp_DIASTOLICadjMEDz", "UKB_460K.bp_SYSTOLICadjMEDz",
                  "PASS_Alzheimer", "PASS_Anorexia", "PASS_Bipolar_Disorder",
                  "PASS_Schizophrenia", "UKB_460K.mental_NEUROTICISM", "UKB_460K.disease_AID_SURE",
                  "PASS_Type_2_Diabetes")

get_sd_annot = function(cell_path, annot_index = 1, base_path, flag=0){
  if(flag == 0){
    sd_annot = rep(0, length(annot_index))
    for(i in 1:length(annot_index)){
      if(file.exists(paste0(cell_path, "/", "sd_annot_", annot_index[i], ".rda"))){
        sd_annot[i] = as.numeric(get(load(paste0(cell_path, "/", "sd_annot_", annot_index[i], ".rda"))))
      }else{
        flag = 1
        break
      }
    }
  }

  if(flag == 1){
    num = rep(0, length(annot_index))
    den = rep(0, length(annot_index))
    ll <- list.files(cell_path, pattern = ".annot.gz")
    ordering = c(1, 10:19, 2, 20:22, 3:9)
    for(m in 1:length(ll)){
      dat <- data.frame(fread(paste0("zcat ", cell_path, "/", ll[m])))
      base <- data.frame(fread(paste0("zcat ", base_path, "/", "baselineLD.", ordering[m], ".annot.gz")))
      pooled_dat <- cbind(dat[,-(1:4)], base[,-(1:4)])
      num = num  + (nrow(pooled_dat)-1) * apply(pooled_dat[,annot_index], 2, var)
      den = den + (nrow(pooled_dat)-1)
      rm(pooled_dat)
    }
    sd_annot = sqrt(num/den)
    for(i in 1:length(annot_index)){
      temp = sd_annot[i]
      save(temp, file = paste0(cell_path, "/", "sd_annot_", annot_index[i], ".rda"))
    }
  }
  return(sd_annot)
}

run_many_tau_analysis = function(annot_cell,
                                 results_cell,
                                 base_path,
                                 annotation,
                                 traits,
                                 index_in_results=NULL,  ### else specify 1:5 for first 5 annotations
                                 base_index = NULL,  ### number of annotations in the baseline
                                 flag = 1){
  base <- data.frame(fread(paste0("zcat ", base_path, "/", "baselineLD.", 22, ".annot.gz")))
  if(is.null(base_index)){
    base_index = ncol(base) - 4
  }
  cell_path = paste0(annot_cell, "/", annotation)
  res = paste0(results_cell, "/", annotation, "/", traits[1], ".sumstats.results")
  tab2 = read.table(res,header=T)
  cat("Number of annotations together with baseline : ", nrow(tab2) , "\n")
  if(is.null(index_in_results)){index_in_results = 1:(nrow(tab2) - base_index)}
  tau_star_table = matrix(0, length(index_in_results), 3)
  annot_names = as.character(tab2$Category[index_in_results])
  sd_annot = get_sd_annot(cell_path, annot_index=index_in_results, base_path, flag = flag)
  for(id in 1:length(index_in_results)){
    sd_annot1=sd_annot[id]
    Mref = 5961159
    df = c()
    for(trait_id in 1:length(traits)){
      result.file=paste0(results_cell, "/", annotation, "/", traits[trait_id], ".sumstats.part_delete")
      new_table=read.table(result.file,header=F)
      sc=c()
      logfile = paste(results_cell, "/", annotation, "/", traits[trait_id],".sumstats.log", sep="")
      log = read.table(logfile,h=F,fill=T)
      h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
      coef1=sd_annot1*Mref/h2g
      for(i in 1:dim(new_table)[1]){
        tau1=as.numeric(new_table[i,index_in_results[id]])
        taus1=tau1*coef1
        #taus1 = tau1
        sc=c(sc,taus1)
        #cat("Block ", i, "\n")
      }
      mean_sc=mean(sc)
      se_sc=sqrt(199**2/200*var(sc))
      df = rbind(df, c(mean_sc,se_sc))
    }
    test_tauj=meta.summaries(df[,1],df[,2],method="random")
    tau=test_tauj$summary
    tau_se=test_tauj$se.summary
    z=tau/tau_se
    cat("Printing results for annotation:", annot_names[id], "\n")
    cat(tau, " ", tau_se, " ", 2*pnorm(-abs(z)), "\n")
    tau_star_table[id, ] = c(tau, tau_se, 2*pnorm(-abs(z)))
  }
  rownames(tau_star_table) = annot_names
  return(tau_star_table)
}


##########################################     brain      ###############################################

out1 = run_many_tau_analysis(annot_cell, results_cell, base_path,
                             annotation = annot_names,
                             traits = all_traits,
                             flag = 1, base_index = 86)
out2 = run_many_tau_analysis(annot_cell, results_cell, base_path,
                             annotation = annot_names,
                             traits = brain_traits,
                             flag = 1, base_index = 106)
out3 = run_many_tau_analysis(annot_cell, results_cell, base_path,
                             annotation = annot_names,
                             traits = c(blood_traits, autoimmune_traits),
                             flag = 0, base_index = NULL)
out4 = run_many_tau_analysis(annot_cell, results_cell, base_path,
                             annotation = annot_names,
                             traits = c(autoimmune_traits),
                             flag = 0, base_index = 86)

total_df = c()
for(name in list.files(results_cell)){
  out3 = run_many_tau_analysis(annot_cell, results_cell, base_path,
                               annotation = name,
                               traits = c(blood_traits, autoimmune_traits),
                               flag = 0, base_index = NULL)
  total_df = rbind(total_df, out3)
  cat("Finished processing files for annot:", name, "\n")
}


ll <- list()
ll[["All"]] = out1
ll[["Brain"]] = out2
ll[["Blood"]] = out3








ll$All[order(ll$All[,3], decreasing = FALSE), ]



out1[which(out1[,3] < 0.1),]

save(ll, file = paste0("/n/groups/price/kushal/LDSC-Average/output/", annot_names, ".rda"))

