# Valentinas Sungaila
# Methods Paper Code


# Graphing an F distribution with df1 =5 and df2=5
# to show how the distribution changes as the noncentrality parameter changes

# Create the vector x
x <- seq(from = 0, to = 9, lenght = 200)

# Evaluate the densities
y_1 <- df(x, 5, 5, ncp = 0)
y_2 <- df(x, 5, 5, ncp = 1)
y_3 <- df(x, 5, 5, ncp = 5)
y_4 <- df(x, 5, 5, ncp = 10)
y_5 <- df(x, 5, 5, ncp = 15)
y_6 <- df(x, 5, 5, ncp = 20)
y_7 <- df(x, 5, 5, ncp = 25)

# Plot the densities
plot(x, y_1, col = 1, type = "l", ylab = '', xlab = '',
     main = "Change in F-Distribution as the Non-Centrality Parameter Changes")
lines(x, y_2, col = 2)
lines(x, y_3, col = 3)
lines(x, y_4, col = 4)
lines(x, y_5, col = 5)
lines(x, y_6, col = 6)
lines(x, y_7, col = 7)

# Add the legend
legend("topright", title = "graph",
       c("df = (5,5) and ncp = 0", "df = (5,5) and ncp = 1", 
         "df = (5,5) and ncp = 5", "df = (5,5) and ncp = 10", 
         "df = (5,5) and ncp = 15", "df = (5,5) and ncp = 20", 
         "df = (5,5) and ncp = 25"), 
       col = c(1, 2, 3, 4, 5, 6, 7), lty = 1)


###-------------------------------------------------------------------------------

# Computing power for F-distribution

# calculate power with 6 different ncp values and changing the number of observations in model

# n = 10
a = 0.05  # type1 error
n1 = 10 # sample size
p1 = 0 # other covariates
d = 1 # number of predictors we are interested in
v = n1 - p1 - d
c = qf(1-a, d, v)
ncp = c(0, 1, 5, 10, 15, 20, 25)
power_ncp.change1 = pf(c, d, v, ncp, lower.tail = FALSE)

# n = 50
n2 = 50 # other covariates
v2 = n2 - p1 - d
c2 = qf(1-a, d, v2)
power_ncp.change2 = pf(c2, d, v2, ncp, lower.tail = FALSE)

# n = 200
n3 = 200 # other covariates
v3 = n3 - p1 - d
c3 = qf(1-a, d, v3)
power_ncp.change3 = pf(c3, d, v3, ncp, lower.tail = FALSE)

# n = 10000
n4 = 10000 # other covariates
v4 = n4 - p1 - d
c4 = qf(1-a, d, v4)
power_ncp.change4 = pf(c4, d, v4, ncp, lower.tail = FALSE)

plot(ncp,power_ncp.change1,
     type = "l",
     main = "Change in power as Non-Centrality Parameter Increases, for Different Number of Observations",
     ylab = "Power",
     xlab = "Non-Centrality Parameter")

lines(ncp, power_ncp.change2, col = 2)
lines(ncp, power_ncp.change3, col = 3)
lines(ncp, power_ncp.change4, col = 4)

# Add the legend
legend("bottomright", title = "graph",
       c("n = 10", "n = 100", 
         "n = 1000", "n = 10000"), 
       col = c(1, 2, 3, 4), lty = 1)

# add the 80% threshold line
abline(h = 0.80, lwd = 2)


###-----------------------------------------------------------------------------------------------------------------------
# graphing # of parameters as expalined by the proportion of the variation of the predictor,
# for power of .8

library(pwr)
u = 1 # # number of predictors we are interested in
# v is the calculated sample size
power = 0.80 # the desired power level

number.of.observations_R2.0.1 = pwr.f2.test(u = u, f2 = 0.1/(1-0.1), sig.level = 0.05, power = power)
number.of.observations_R2.0.2 = pwr.f2.test(u = u, f2 = 0.2/(1-0.2), sig.level = 0.05, power = power)
number.of.observations_R2.0.3 = pwr.f2.test(u = u, f2 = 0.3/(1-0.3), sig.level = 0.05, power = power)
number.of.observations_R2.0.4 = pwr.f2.test(u = u, f2 = 0.4/(1-0.4), sig.level = 0.05, power = power)
number.of.observations_R2.0.5 = pwr.f2.test(u = u, f2 = 0.5/(1-0.5), sig.level = 0.05, power = power)
number.of.observations_R2.0.6 = pwr.f2.test(u = u, f2 = 0.6/(1-0.6), sig.level = 0.05, power = power)
number.of.observations_R2.0.7 = pwr.f2.test(u = u, f2 = 0.7/(1-0.7), sig.level = 0.05, power = power)
number.of.observations_R2.0.8 = pwr.f2.test(u = u, f2 = 0.8/(1-0.8), sig.level = 0.05, power = power)
number.of.observations_R2.0.9 = pwr.f2.test(u = u, f2 = 0.9/(1-0.9), sig.level = 0.05, power = power)

number_of_observations_per_R2 = rbind("R2 = 0.01" = number.of.observations_R2.0.1$v,
                                      "R2 = 0.02" = number.of.observations_R2.0.2$v,
                                      "R2 = 0.03" = number.of.observations_R2.0.3$v,
                                      "R2 = 0.04" = number.of.observations_R2.0.4$v,
                                      "R2 = 0.05" = number.of.observations_R2.0.5$v,
                                      "R2 = 0.06" = number.of.observations_R2.0.6$v,
                                      "R2 = 0.07" = number.of.observations_R2.0.7$v,
                                      "R2 = 0.08" = number.of.observations_R2.0.8$v,
                                      "R2 = 0.09" = number.of.observations_R2.0.9$v)

number_of_observations_per_R2= cbind(number_of_observations_per_R2,
                                    c(0.01,
                                     0.02,
                                     0.03,
                                     0.04,
                                     0.05,
                                     0.06,
                                     0.07,
                                     0.08,
                                     0.09))


number_of_observations_per_R2_1 = as.data.frame(number_of_observations_per_R2)

#####
u = 10 # # number of predictors we are interested in
# v is the calculated sample size
power = 0.80 # the desired power level

number.of.observations_R2.0.1 = pwr.f2.test(u = u, f2 = 0.1/(1-0.1), sig.level = 0.05, power = power)
number.of.observations_R2.0.2 = pwr.f2.test(u = u, f2 = 0.2/(1-0.2), sig.level = 0.05, power = power)
number.of.observations_R2.0.3 = pwr.f2.test(u = u, f2 = 0.3/(1-0.3), sig.level = 0.05, power = power)
number.of.observations_R2.0.4 = pwr.f2.test(u = u, f2 = 0.4/(1-0.4), sig.level = 0.05, power = power)
number.of.observations_R2.0.5 = pwr.f2.test(u = u, f2 = 0.5/(1-0.5), sig.level = 0.05, power = power)
number.of.observations_R2.0.6 = pwr.f2.test(u = u, f2 = 0.6/(1-0.6), sig.level = 0.05, power = power)
number.of.observations_R2.0.7 = pwr.f2.test(u = u, f2 = 0.7/(1-0.7), sig.level = 0.05, power = power)
number.of.observations_R2.0.8 = pwr.f2.test(u = u, f2 = 0.8/(1-0.8), sig.level = 0.05, power = power)
number.of.observations_R2.0.9 = pwr.f2.test(u = u, f2 = 0.9/(1-0.9), sig.level = 0.05, power = power)

number_of_observations_per_R2 = rbind("R2 = 0.01" = number.of.observations_R2.0.1$v,
                                      "R2 = 0.02" = number.of.observations_R2.0.2$v,
                                      "R2 = 0.03" = number.of.observations_R2.0.3$v,
                                      "R2 = 0.04" = number.of.observations_R2.0.4$v,
                                      "R2 = 0.05" = number.of.observations_R2.0.5$v,
                                      "R2 = 0.06" = number.of.observations_R2.0.6$v,
                                      "R2 = 0.07" = number.of.observations_R2.0.7$v,
                                      "R2 = 0.08" = number.of.observations_R2.0.8$v,
                                      "R2 = 0.09" = number.of.observations_R2.0.9$v)

number_of_observations_per_R2= cbind(number_of_observations_per_R2,
                                     c(0.01,
                                       0.02,
                                       0.03,
                                       0.04,
                                       0.05,
                                       0.06,
                                       0.07,
                                       0.08,
                                       0.09))


number_of_observations_per_R2_2 = as.data.frame(number_of_observations_per_R2)


#####
u = 100 # # number of predictors we are interested in
# v is the calculated sample size
power = 0.80 # the desired power level

number.of.observations_R2.0.1 = pwr.f2.test(u = u, f2 = 0.1/(1-0.1), sig.level = 0.05, power = power)
number.of.observations_R2.0.2 = pwr.f2.test(u = u, f2 = 0.2/(1-0.2), sig.level = 0.05, power = power)
number.of.observations_R2.0.3 = pwr.f2.test(u = u, f2 = 0.3/(1-0.3), sig.level = 0.05, power = power)
number.of.observations_R2.0.4 = pwr.f2.test(u = u, f2 = 0.4/(1-0.4), sig.level = 0.05, power = power)
number.of.observations_R2.0.5 = pwr.f2.test(u = u, f2 = 0.5/(1-0.5), sig.level = 0.05, power = power)
number.of.observations_R2.0.6 = pwr.f2.test(u = u, f2 = 0.6/(1-0.6), sig.level = 0.05, power = power)
number.of.observations_R2.0.7 = pwr.f2.test(u = u, f2 = 0.7/(1-0.7), sig.level = 0.05, power = power)
number.of.observations_R2.0.8 = pwr.f2.test(u = u, f2 = 0.8/(1-0.8), sig.level = 0.05, power = power)
number.of.observations_R2.0.9 = pwr.f2.test(u = u, f2 = 0.9/(1-0.9), sig.level = 0.05, power = power)

number_of_observations_per_R2 = rbind("R2 = 0.01" = number.of.observations_R2.0.1$v,
                                      "R2 = 0.02" = number.of.observations_R2.0.2$v,
                                      "R2 = 0.03" = number.of.observations_R2.0.3$v,
                                      "R2 = 0.04" = number.of.observations_R2.0.4$v,
                                      "R2 = 0.05" = number.of.observations_R2.0.5$v,
                                      "R2 = 0.06" = number.of.observations_R2.0.6$v,
                                      "R2 = 0.07" = number.of.observations_R2.0.7$v,
                                      "R2 = 0.08" = number.of.observations_R2.0.8$v,
                                      "R2 = 0.09" = number.of.observations_R2.0.9$v)

number_of_observations_per_R2= cbind(number_of_observations_per_R2,
                                     c(0.01,
                                       0.02,
                                       0.03,
                                       0.04,
                                       0.05,
                                       0.06,
                                       0.07,
                                       0.08,
                                       0.09))


number_of_observations_per_R2_3 = as.data.frame(number_of_observations_per_R2)



#####
u = 1000 # # number of predictors we are interested in
# v is the calculated sample size
power = 0.80 # the desired power level

number.of.observations_R2.0.1 = pwr.f2.test(u = u, f2 = 0.1/(1-0.1), sig.level = 0.05, power = power)
number.of.observations_R2.0.2 = pwr.f2.test(u = u, f2 = 0.2/(1-0.2), sig.level = 0.05, power = power)
number.of.observations_R2.0.3 = pwr.f2.test(u = u, f2 = 0.3/(1-0.3), sig.level = 0.05, power = power)
number.of.observations_R2.0.4 = pwr.f2.test(u = u, f2 = 0.4/(1-0.4), sig.level = 0.05, power = power)
number.of.observations_R2.0.5 = pwr.f2.test(u = u, f2 = 0.5/(1-0.5), sig.level = 0.05, power = power)
number.of.observations_R2.0.6 = pwr.f2.test(u = u, f2 = 0.6/(1-0.6), sig.level = 0.05, power = power)
number.of.observations_R2.0.7 = pwr.f2.test(u = u, f2 = 0.7/(1-0.7), sig.level = 0.05, power = power)
number.of.observations_R2.0.8 = pwr.f2.test(u = u, f2 = 0.8/(1-0.8), sig.level = 0.05, power = power)
number.of.observations_R2.0.9 = pwr.f2.test(u = u, f2 = 0.9/(1-0.9), sig.level = 0.05, power = power)

number_of_observations_per_R2 = rbind("R2 = 0.01" = number.of.observations_R2.0.1$v,
                                      "R2 = 0.02" = number.of.observations_R2.0.2$v,
                                      "R2 = 0.03" = number.of.observations_R2.0.3$v,
                                      "R2 = 0.04" = number.of.observations_R2.0.4$v,
                                      "R2 = 0.05" = number.of.observations_R2.0.5$v,
                                      "R2 = 0.06" = number.of.observations_R2.0.6$v,
                                      "R2 = 0.07" = number.of.observations_R2.0.7$v,
                                      "R2 = 0.08" = number.of.observations_R2.0.8$v,
                                      "R2 = 0.09" = number.of.observations_R2.0.9$v)

number_of_observations_per_R2= cbind(number_of_observations_per_R2,
                                     c(0.01,
                                       0.02,
                                       0.03,
                                       0.04,
                                       0.05,
                                       0.06,
                                       0.07,
                                       0.08,
                                       0.09))


number_of_observations_per_R2_4 = as.data.frame(number_of_observations_per_R2)

#####

plot(number_of_observations_per_R2_1$V2, number_of_observations_per_R2_1$V1,
     type = "l",
     xlab = "R^2",
     ylab = "Number of Observations",
     main = "Change in the Number of Observations as R-squared Increases for Power = 0.8")


lines(number_of_observations_per_R2_1$V2, number_of_observations_per_R2_2$V1, col = 2)
lines(number_of_observations_per_R2_1$V2, number_of_observations_per_R2_3$V1, col = 3)
lines(number_of_observations_per_R2_1$V2, number_of_observations_per_R2_4$V1, col = 4)

# Add the legend
legend("topright", title = "graph",
       c("# of predictors = 1", "# of predictors = 10", 
         "# of predictors = 100", "# of predictors = 1000"), 
       col = c(1, 2, 3, 4), lty = 1)


######_________________________________________________________________________________________________

######### Libraries #########
library(data.table)
library(nnet) 
library(rms)
library(lmtest)
library(pwr)
library(factoextra)



######### Data #########
dir_path <- "C:/Users/sunga/Desktop/consulting/" #WALDO, change me
dat <- as.data.frame(fread(paste0(dir_path,"subphenAndGenoVarsToUse.csv")))

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












