# Genetic variants and vitiligo disease characteristics
# Authors: Sandra Robles Munoz & Valentinas Sungaila
# Spring 2020, MATH 6330

######### Libraries #########
library(data.table)
library(nnet) 
library(rms)
library(lmtest)
library(pwr)
library(factoextra)



######### Data #########
dir_path <- "~/Grad School/MATH-6330_S20/Project/" #WALDO, change me
dat <- as.data.frame(fread(paste0(dir_path,"Data_6330_Project.csv")))

######### Data Hygiene #########
colnames(dat) <- gsub("bcsnp.vitiligo_assoc","vit_ass",colnames(dat))
colnames(dat) <- gsub("bcsnp.vitiligo","vit_ass",colnames(dat))
colnames(dat) <- gsub("bcsnp.patient","patient",colnames(dat))

dat$vit_ass.rheumatoid_arthritis <- as.integer(ifelse(dat$vit_ass.rheumatoid_arthritis=='yes',1,0))
dat$vit_ass.systemic_lupus_erythematosus <- as.integer(ifelse(dat$vit_ass.systemic_lupus_erythematosus=='yes',1,0))
dat$vit_ass.pernicious_anemia <- as.integer(ifelse(dat$vit_ass.pernicious_anemia=='yes',1,0))
dat$vit_ass_onset.acrofacial <- as.integer(ifelse(dat$vit_ass_onset.acrofacial=='yes',1,0))
dat$vit_ass.halo_nevi <- as.integer(ifelse(dat$vit_ass.halo_nevi=='yes',1,0))
dat$vit_ass_onset.koebner_phenomenon <- as.integer(ifelse(dat$vit_ass_onset.koebner_phenomenon=='yes',1,0))

genotype_cols <- grep('^RS',colnames(dat),value=T)
covariates_cols <- grep('^GWAS',colnames(dat),value=T)
overlay_cols <- grep('patient.|vit_ass',colnames(dat),value=T)
dependent_cols <- grep('SKIN|DURATION',colnames(dat),value=T)
score_cols <- grep("^score",colnames(dat),value = T)

dat$SKIN <- factor(dat$SKIN,ordered=TRUE,levels=c('up to 25%','26-50%','51-75%','76-100%'))

######### PRE-PROCESSING: Picking top Ancestry Principal Components #########

# All gwas are populated otherwise we'd have to filter down to full observations
fwer_value <- 0.10/length(covariates_cols) # Using family-wise error adjustment (alpha/number of tests)
ageonset_keep <- c()
skin_keep <- c()

# Step 1 - Duration 
mean_model_age <- lm(patient.vitiligo_age_of_onset ~ 1, data = dat)
for(count in covariates_cols){
  temp_model <- lm(paste0("patient.vitiligo_age_of_onset ~ ",count),data = dat)
  fwer_test <- anova(mean_model_age,temp_model)$`Pr(>F)`[2] < fwer_value
  if(fwer_test){ageonset_keep <- c(ageonset_keep,count)}
}

# Step 2 - SKIN

# No need to create a mean model, the 'orm' function natively does a loglikelihood test
for(count in covariates_cols){
  temp_model <- orm(as.formula(paste0("SKIN ~ ",count)),data = dat)
  fwer_test <- as.numeric(temp_model$stats[7]) < fwer_value
  print(as.numeric(temp_model$stats[7]))
  if(fwer_test){skin_keep <- c(skin_keep,count)}
  # Note: no PCs kept using ordinal model, GWAS1.EV1 kept if nominal.However, ordinality is respected, so orm() is used instead of multinom()
  # Note 2: POLR fails on GWAS1.V1 due to 'initial value in 'vmmin' is not finite', google search indicates that using the package RMS bypasses this problem
  
}


######### QUESTION 1 - Duration #########
# Filters down to complete observations
age_model_list <- list()
age_df <- data.frame(genetic_variant=character(),
                     Coefficient=double(),
                    var_prob=double(),
                    std_err=double(),
                    Rsqd=double(),
                    Rsqrd_adj=double(),
                    pval_anova=double(),stringsAsFactors = F)
if(length(ageonset_keep>0)){ #PC covariates will be applied
 for(count in genotype_cols){
   age_model_list[[count]] <- lm(paste0("patient.vitiligo_age_of_onset ~ ",count," + ",paste0(ageonset_keep,collapse = " + ")),data = dat[!is.na(colnames(dat) %in% c(count))])
   tmp_summary <- summary(age_model_list[[count]])
   tmp_coeff <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Estimate"]
   tmp_stderr <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Std. Error"]
   tmp_valprob <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Pr(>|t|)"]
   tmp_rsqrd <- tmp_summary$r.squared
   tmp_rsqrdadj <- tmp_summary$adj.r.squared
   tmp_fstatprob <- as.double(pf(tmp_summary$fstatistic[1],tmp_summary$fstatistic[2],tmp_summary$fstatistic[3],lower.tail=FALSE))
   age_df[nrow(age_df)+1,]<- c(count,tmp_coeff,
                               tmp_valprob,
                               tmp_stderr,
                               tmp_rsqrd,
                               tmp_rsqrdadj,
                               tmp_fstatprob)
 }
}else{
  for(count in genotype_cols){
    age_model_list[[count]] <- lm(paste0("patient.vitiligo_age_of_onset ~ ",count),data = dat[!is.na(colnames(dat) %in% c(count))])
    tmp_summary <- summary(age_model_list[[count]])
    tmp_coeff <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Estimate"]
    tmp_stderr <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Std. Error"]
    tmp_valprob <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Pr(>|t|)"]
    tmp_rsqrd <- as.double(tmp_summary$r.squared)
    tmp_rsqrdadj <- as.double(tmp_summary$adj.r.squared)
    tmp_fstatprob <- as.double(pf(tmp_summary$fstatistic[1],tmp_summary$fstatistic[2],tmp_summary$fstatistic[3],lower.tail=FALSE))
    age_df[nrow(age_df)+1,]<- c(count,tmp_coeff,
                                tmp_valprob,
                                tmp_stderr,
                                tmp_rsqrd,
                                tmp_rsqrdadj,
                                tmp_fstatprob)
  }
  
}
age_df[!colnames(age_df) %in% c("genetic_variant")] <- lapply(age_df[!colnames(age_df) %in% c("genetic_variant")], as.numeric)

######### QUESTION 1 - Skin surface #########
# The model can be built with skin treated as ordinal or nominal. There's often a lot of information to be gained from leaving alone the ordinal nature of the outcome

skin_model_list <- list()
skin_df <- data.frame(genetic_variant=character(),
                      likelihoodration=double(),
                      df=integer(),
                      likelihood_pval=double(),
                      chisqrd_val=double(),
                      chisqrd_pval=double(),
                      stringsAsFactors = F)
if(length(skin_keep>0)){ #PC covariates will be applied
  for(count in genotype_cols){
    skin_model_list[[count]] <- orm(as.formula(paste0("SKIN ~ ",count," + ",paste0(skin_keep,collapse = " + "))),data=dat[!is.na(colnames(dat) %in% c(count))])
    tmp_summary <- skin_model_list[[count]]$stats
    tmp_likelihoodratio = as.double(tmp_summary[5])
    tmp_df = as.integer(tmp_summary[6])
    tmp_LR_pval = as.double(tmp_summary[7])
    tmp_chisqrd = as.double(tmp_summary[8])
    tmp_chisqrd_pval = as.double(tmp_summary[9])
    #pchisq(as.double(tmp_summary[8]),df = 2,lower.tail = F)
    skin_df[nrow(skin_df)+1,]<- c(count,tmp_likelihoodratio,
                                  tmp_df,
                                  tmp_LR_pval,
                                  tmp_chisqrd,tmp_chisqrd_pval)
  }
  
}else{
  for(count in genotype_cols){
    skin_model_list[[count]] <- orm(as.formula(paste0("SKIN ~ ",count)),data=dat[!is.na(colnames(dat) %in% c(count))])
    tmp_summary <- skin_model_list[[count]]$stats
    tmp_likelihoodratio = as.double(tmp_summary[5])
    tmp_df = as.integer(tmp_summary[6])
    tmp_LR_pval = as.double(tmp_summary[7])
    tmp_chisqrd = as.double(tmp_summary[8])
    tmp_chisqrd_pval = as.double(tmp_summary[9])
    #pchisq(as.double(tmp_summary[8]),df = 2,lower.tail = F)
    skin_df[nrow(skin_df)+1,]<- c(count,tmp_likelihoodratio,
                                  tmp_df,
                                  tmp_LR_pval,
                                  tmp_chisqrd,tmp_chisqrd_pval)
  }
  
}
skin_df[!colnames(skin_df) %in% c("genetic_variant")] <- lapply(skin_df[!colnames(skin_df) %in% c("genetic_variant")], as.numeric)





######### QUESTION 2 - Duration #########
# Q2: What effect sizes for disease characteristics are detectable with >= 80% power 
        # without adjustment for multiple testing, 
        # with adjust for the number of tests in Q1, 
        # and with adjustment genome-wide?

n_ageonset <- nrow(dat[!is.na(dat$patient.vitiligo_age_of_onset),])
df_num = 2-1 # Predictors minus intercept. We don't have any covariates to account for, otherwise that'd go here
effect_unadj <- pwr.f2.test(u = df_num, v=n_ageonset-1-1, sig.level = 0.05, power = 0.8)
effect_bonferroniadj <- pwr.f2.test(u = df_num, v=n_ageonset-1-1, sig.level = 0.05/46, power = 0.8)
effect_genomewide <- pwr.f2.test(u = df_num, v=n_ageonset-1-1, sig.level = 5e-8, power = 0.8)

# Effect (R^2) needed for 80% power with a 2190 sample size:
eff_unadj <- effect_unadj$f2/(1+effect_unadj$f2)  #0.003576067
eff_adj <- effect_bonferroniadj$f2/(1+effect_bonferroniadj$f2) #0.007686651
eff_genadj <- effect_genomewide$f2/(1+effect_genomewide$f2) #0.01788534

age_df$effect_largeenough_unadj <- ifelse(age_df$Rsqd>=eff_unadj,'Sufficient effect','Insufficient effect')
age_df$effect_largeenough_bonfadj <- ifelse(age_df$Rsqd>=eff_adj,'Sufficient effect','Insufficient effect')
age_df$effect_largeenough_genomeadj<- ifelse(age_df$Rsqd>=eff_genadj,'Sufficient effect','Insufficient effect')
# Only RS114448410 has a big enough effect size for 80% power at the unadjusted and bonferroni adjustment level. None at genone-wide


######### QUESTION 2 - Skin surface #########
#NOTE: Power is not readily available for Multinomial Logistic Regression and its assumed distribution
# So this step will requite recreating the 46 regressions with SKIN recoded as numerical, and then run as linear models
#dat_skin <- dat[!is.na(dat$SKIN),]
dat$SKIN_numeric <- ifelse(dat$SKIN=='up to 25%',1,
      ifelse(dat$SKIN=='26-50%',2,
             ifelse(dat$SKIN=='51-75%',3,
                    ifelse(dat$SKIN=='76-100%',4,dat$SKIN
       )))) 

# First figure out any covariates
skin_numerical_keep <- c()
mean_model_ageNumerical <- lm(SKIN_numeric ~ 1, data = dat)
for(count in covariates_cols){
  temp_model <- lm(paste0("SKIN_numeric ~ ",count),data = dat)
  fwer_test <- anova(mean_model_ageNumerical,temp_model)$`Pr(>F)`[2] < fwer_value
  if(fwer_test){skin_numerical_keep <- c(skin_numerical_keep,count)}
}
# Now we construct the linear models to get the R^2 values
skinnumerical_model_list <- list()
skinnumerical_df <- data.frame(genetic_variant=character(),
                               Coefficient=double(),
                               var_prob=double(),
                               std_err=double(),
                               Rsqd=double(),
                               Rsqrd_adj=double(),
                               pval_anova=double(),stringsAsFactors = F)
if(length(ageonset_keep>0)){ #PC covariates will be applied
  for(count in genotype_cols){
    skinnumerical_model_list[[count]] <- lm(paste0("SKIN_numeric ~ ",count," + ",paste0(ageonset_keep,collapse = " + ")),data = dat[!is.na(colnames(dat) %in% c(count))])
    tmp_summary <- summary(skinnumerical_model_list[[count]])
    tmp_coeff <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Estimate"]
    tmp_stderr <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Std. Error"]
    tmp_valprob <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Pr(>|t|)"]
    tmp_rsqrd <- tmp_summary$r.squared
    tmp_rsqrdadj <- tmp_summary$adj.r.squared
    tmp_fstatprob <- as.double(pf(tmp_summary$fstatistic[1],tmp_summary$fstatistic[2],tmp_summary$fstatistic[3],lower.tail=FALSE))
    skinnumerical_df[nrow(skinnumerical_df)+1,]<- c(count,tmp_coeff,
                                                    tmp_valprob,
                                                    tmp_stderr,
                                                    tmp_rsqrd,
                                                    tmp_rsqrdadj,
                                                    tmp_fstatprob)
  }
}else{
  for(count in genotype_cols){
    skinnumerical_model_list[[count]] <- lm(paste0("SKIN_numeric ~ ",count),data = dat[!is.na(colnames(dat) %in% c(count))])
    tmp_summary <- summary(skinnumerical_model_list[[count]])
    tmp_coeff <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Estimate"]
    tmp_stderr <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Std. Error"]
    tmp_valprob <- tmp_summary$coefficients[row.names(tmp_summary$coefficients)==count,"Pr(>|t|)"]
    tmp_rsqrd <- as.double(tmp_summary$r.squared)
    tmp_rsqrdadj <- as.double(tmp_summary$adj.r.squared)
    tmp_fstatprob <- as.double(pf(tmp_summary$fstatistic[1],tmp_summary$fstatistic[2],tmp_summary$fstatistic[3],lower.tail=FALSE))
    skinnumerical_df[nrow(skinnumerical_df)+1,]<- c(count,tmp_coeff,
                                                    tmp_valprob,
                                                    tmp_stderr,
                                                    tmp_rsqrd,
                                                    tmp_rsqrdadj,
                                                    tmp_fstatprob)
  }
  
}
skinnumerical_df[!colnames(skinnumerical_df) %in% c("genetic_variant")] <- lapply(skinnumerical_df[!colnames(skinnumerical_df) %in% c("genetic_variant")], as.numeric)


# And now the Power calculations:
n_skin <- nrow(dat[!is.na(dat$SKIN_numeric),])
df_num = 2-1 # Predictors minus intercept. We don't have any covariates to account for, otherwise that'd go here
effect_unadj_skin <- pwr.f2.test(u = df_num, v=n_skin-1-1, sig.level = 0.05, power = 0.8)
effect_bonferroniadj_skin <- pwr.f2.test(u = df_num, v=n_skin-1-1, sig.level = 0.05/46, power = 0.8)
effect_genomewide_skin <- pwr.f2.test(u = df_num, v=n_skin-1-1, sig.level = 5e-8, power = 0.8)

# Effect (R^2) needed for 80% power with a 1988 sample size:
eff_unadj_skin <- effect_unadj_skin$f2/(1+effect_unadj_skin$f2)  # 0.00395165
eff_adj_skin <- effect_bonferroniadj_skin$f2/(1+effect_bonferroniadj_skin$f2) # 0.008435601
eff_genadj_skin <- effect_genomewide_skin$f2/(1+effect_genomewide_skin$f2) # 0.02008278

skinnumerical_df$effect_largeenough_unadj <- ifelse(skinnumerical_df$Rsqd>=eff_unadj_skin,'Sufficient effect','Insufficient effect')
skinnumerical_df$effect_largeenough_bonfadj <- ifelse(skinnumerical_df$Rsqd>=eff_adj_skin,'Sufficient effect','Insufficient effect')
skinnumerical_df$effect_largeenough_genomeadj<- ifelse(skinnumerical_df$Rsqd>=eff_genadj_skin,'Sufficient effect','Insufficient effect')
# Only RS10986311 has a big enough effect size for 80% power at the unadjusted level. None at bonferroni or genone-wide level.



######### QUESTION 3 - Structure of cases based on disease characteristics? #########
# structure of cases can be derived using clustering algorithms, however,  they cannot handle categorical data or missing data
# For categorical: clustering depends on distances (dissimilarities) which is almost ways done with eucledian distances
# For missing: same reason of distance calculations. There are imputation techniques available

cluster_data <- dat[,c(score_cols,'SKIN_numeric','patient.vitiligo_age_of_onset')]
cluster_data <- na.omit(cluster_data)
cluster_data <- scale(cluster_data)

clustersss <-kmeans(x=as.matrix(cluster_data),centers = 3,nstart = 150, iter.max=30)
table(clustersss$cluster)
fviz_cluster(clustersss, geom = "point", data = cluster_data) + ggtitle("First 2 principal components by cluster")



######### PRESENTATION/PAPER CODE ######### 
dat$temp_ageonset = ifelse(between(dat$patient.vitiligo_age_of_onset,0,10),'0-10',
                           ifelse(between(dat$patient.vitiligo_age_of_onset,11,20),'11-20',
                                  ifelse(between(dat$patient.vitiligo_age_of_onset,21,30),'21-30',
                                         ifelse(between(dat$patient.vitiligo_age_of_onset,31,40),'31-40',
                                                ifelse(between(dat$patient.vitiligo_age_of_onset,41,50),'41-50',
                                                       ifelse(between(dat$patient.vitiligo_age_of_onset,51,60),'51-60',
                                                              ifelse(between(dat$patient.vitiligo_age_of_onset,61,70),'61-70',
                                                                     ifelse(between(dat$patient.vitiligo_age_of_onset,71,80),'71-80',
                                                                            ifelse(dat$patient.vitiligo_age_of_onset>80,'80+',dat$patient.vitiligo_age_of_onset)))))))))
table(dat$SKIN,dat$temp_ageonset)
 