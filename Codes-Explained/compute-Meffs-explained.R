#COPD Gene LC-MS Data
#Simulation scenario: sub-sample size nsub
#number of metabolomic features p=p
#Outcome: log-FEV-ratio -> logarithm of r/(1-r), where r = FEV1/FVC
#FEV1 = Forced Expiratory Volume at 1 second
#FVC = Forced Vital Capacity

#Install package "limma"
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("limma")

library(limma) #load the Bioconductor package "limma" in the current session
library(devtools) #load the R-package "devtools" in the current session
#install the GitHub sub-directory "MWSL" corresponding to the article Peluso et al. (2021), BMC Bioinformatics volume 22, Article number: 67 
install_github("AlinaPeluso/PhenoMeNal", subdir="MWSL") 
library(MWSL) #load the R-package "MWSL" in the current session

rm(list=ls()) #remove all objects from the current session

#Load the final, cleaned metabolomics data with all subjects and p metabolites stored within the RData object "final_processed_dataframes_subjects_metabolites.RData"
load("final_processed_dataframes_subjects_metabolites.RData")


#Select REP-many random subsets, each of size nsub, from the total size of n=1054 patients
nsub = n #sub-sample size
REP = R #no of sub-samples
sampSel=matrix(NA, REP, nsub) #initial a matrix of size REP-by-nsub, where each element = NA
#initial three matrices each of size REP-by-2, where each element = NA
genderSel = smokeSel = ccenterSel = matrix(NA, REP, 2) 

#Check the distribution of the ratio r=FEV1/FVC 
hist(fevrat.final2) #plot the histogram
summary(fevrat.final2) #show the 5-number summary of the continuous variable r
#Min.    1st Qu.  Median    Mean    3rd Qu.    Max. 
#0.2000  0.5700   0.7000    0.6584  0.7700     0.9400

#Store the sample indices of the randomly selected nsub sample-units and the frequency distributions of the categorical covariates, for each of the REP repetitions
for(jj in 1:REP){
  print(jj)
  set.seed(jj)
  ind = sample(length(fevrat.final2), nsub, replace=F)
  
  sampSel[jj,] = ind
  
  #store the frequency distribution of gender (females indexed by 0 and males indexed by 1)
  genderSel[jj,] = table(conf_logFEVrat.final2$GENDER[ind]) 
  
  #store the frequency distribution of smokers (former indexed by 0 and current indexed by 1)
  smokeSel[jj,] = table(conf_logFEVrat.final2$smoking_status[ind])
  
  #store the frequency distribution of sample collection center (NJC indexed by 0 and UIA indexed by 1)
  ccenterSel[jj,] = table(conf_logFEVrat.final2$ccenter[ind])
}

#Save the sample indices for the randomly selected nsub units, for each of the REP repetitions
write.table(sampSel, file = "sampSel_nsub-p_repR.csv", row.names = F, col.names = F)

#Check if the selected nsub sample-units in REP=R repetitions represent 
#the three categorical covariates fairly in their binary categories
par(mfrow=c(1,3))

#minimum of all the REP proportions of females, former smokers, and NJC sample collection center
mini = min(genderSel[,1]/nsub, smokeSel[,1]/nsub, ccenterSel[,1]/nsub)
#maximum of all the REP proportions for females, former smokers, and NJC sample collection center
maxi = max(smokeSel[,1]/nsub, ccenterSel[,1]/nsub, genderSel[,1]/nsub)

#Boxplot of the proportions of females
boxplot(genderSel[,1]/nsub, col=5, 
        ylim = c(mini,maxi),
        ylab="Proportion",
        main="Boxplot of\nProportions of females")

#Boxplot of the proportions of former smokers
boxplot(smokeSel[,1]/nsub, col=3,
        ylim = c(mini,maxi),
        ylab="Proportion",
        main="Boxplot of\nProportions of former smokers")

#Boxplot of the proportions of NJC sample-collection center
boxplot(ccenterSel[,1]/nsub, col=4, 
        ylim = c(mini,maxi),
        ylab="Proportion",
        main="Boxplot of\nProportions of NJC patients")

dev.off()

#Plot all the three boxplots displayed side-by-side within one image
#Concatenate the three proportions in three columns
cat_3 = data.frame(genderSel[,1]/nsub, smokeSel[,1]/nsub, ccenterSel[,1]/nsub)
dim(cat_3) #Display the dimension of the concatenated data-frame
colnames(cat_3) = c("Females", "Former-smokers", "NJC-patients") #Name the columns

#Save the plot as a PNG
png("catCovsProportions_logFEVrat_nsub-p_repR.png",
    width = 12, height = 8, units = "in", res = 300)
boxplot(genderSel[,1]/nsub, smokeSel[,1]/nsub, ccenterSel[,1]/nsub,
        ylim = c(mini,maxi),
        ylab="Proportions",
        xlab = "Categorical covariates",
        names = colnames(cat_3),
        col = rainbow(3),
        main="Proportions of Females, Former Smokers, and NJC-Patients")
dev.off()

#Create lists of corresponding sub-sampled outcome and features (log of metabolite abundances), adjusted-for-covariates features, 
#features into super-pathway groups, confounders 
#Define lists for outcome and covariate matrix, each of length REP
outcome_logFEVrat_sub = conf_sub_scaled = vector(mode = "list", length = REP)

#Define lists for the adjusted-for-covariates metabolites, for each of the 8 sub-groups, each of length REP
logmet_class1_sub_res = logmet_class2_sub_res = logmet_class3_sub_res = logmet_class4_sub_res = 
  logmet_class5_sub_res = logmet_class6_sub_res = logmet_class7_sub_res = logmet_class8_sub_res = 
  vector(mode = "list", length = REP)

for(k in 1:REP){
  print(k)
  
  #collect the corresponding outcome values using the previously stored sample indices
  outcome_logFEVrat_sub[[k]] = fevrat.final2[sampSel[k,]] 
  #collect the corresponding covariate values using the previously stored sample indices
  conf_sub = conf_logFEVrat.final2[sampSel[k,],]
  
  #Center/scale the continuous covariates
  conf_sub_scaled[[k]] = conf_sub
  #Scale the three continuous covariates (BMI, age, ATS pack-years) marginally to zero-mean and unit-variance
  conf_sub_scaled[[k]][, c(1,2,4)] = scale(conf_sub[, c(1,2,4)])
  
  mat = matrix(NA, nsub, ncol(mets.p761.final2)) #initiate matrix of size nsub-by-p 
  
  #store the residuals of metabolites adjusted for the covariates
  for(l in 1:ncol(mets.p761.final2)){
    lin.reg = lm(log(mets.p761.final2[sampSel[k,],l]) ~., data = conf_sub_scaled[[k]])
    mat[,l] = lin.reg$residuals
  }
  
  #save the adjusted-for-covariates metabolite abundances
  filename = paste(c("logfeat_sub_res_nsub-p_", k, ".csv"), collapse = "")
  write.table(mat, filename, row.names = F, col.names = F, sep = ",")
  
  #Create the adjusted-for-covariates metabolite abundances for each of the 8 sub-groups
  logmet_class1_sub_res[[k]] = mat[, 1:365]   #Lipids = 365
  logmet_class2_sub_res[[k]] = mat[, 366:464] #Xenobiotics = 99
  logmet_class3_sub_res[[k]] = mat[, 465:645] #Amino Acid = 181
  logmet_class4_sub_res[[k]] = mat[, 646:669] #Carbohydrate = 24
  logmet_class5_sub_res[[k]] = mat[, 670:694] #Cofactors & Vitamins = 25
  logmet_class6_sub_res[[k]] = mat[, 695:726] #Nucleotide = 32
  logmet_class7_sub_res[[k]] = mat[, 727:736] #Energy = 10
  logmet_class8_sub_res[[k]] = mat[, 737:761] #Peptide = 25
  
  filename1 = paste(c("logfeat_sub_res_nsub-p-class1_", k, ".csv"), collapse = "")
  filename2 = paste(c("logfeat_sub_res_nsub-p-class2_", k, ".csv"), collapse = "")
  filename3 = paste(c("logfeat_sub_res_nsub-p-class3_", k, ".csv"), collapse = "")
  filename4 = paste(c("logfeat_sub_res_nsub-p-class4_", k, ".csv"), collapse = "")
  filename5 = paste(c("logfeat_sub_res_nsub-p-class5_", k, ".csv"), collapse = "")
  filename6 = paste(c("logfeat_sub_res_nsub-p-class6_", k, ".csv"), collapse = "")
  filename7 = paste(c("logfeat_sub_res_nsub-p-class7_", k, ".csv"), collapse = "")
  filename8 = paste(c("logfeat_sub_res_nsub-p-class8_", k, ".csv"), collapse = "")
  
  #Save the lists of 8 metabolite sub-groups, metabolites adjusted for covariates for each sub-group 
  write.table(logmet_class1_sub_res[[k]], filename1, row.names = F, col.names = F, sep = ",")
  write.table(logmet_class2_sub_res[[k]], filename2, row.names = F, col.names = F, sep = ",")
  write.table(logmet_class3_sub_res[[k]], filename3, row.names = F, col.names = F, sep = ",")
  write.table(logmet_class4_sub_res[[k]], filename4, row.names = F, col.names = F, sep = ",")
  write.table(logmet_class5_sub_res[[k]], filename5, row.names = F, col.names = F, sep = ",")
  write.table(logmet_class6_sub_res[[k]], filename6, row.names = F, col.names = F, sep = ",")
  write.table(logmet_class7_sub_res[[k]], filename7, row.names = F, col.names = F, sep = ",")
  write.table(logmet_class8_sub_res[[k]], filename8, row.names = F, col.names = F, sep = ",")
}

#Save the outcome subsets, covariate-adjusted-metabolite abundances (for all p=761 and 8 subgroups)
save(outcome_logFEVrat_sub, conf_sub_scaled, sampSel, mets.p761.final2,
     file = "nsub-p_REP_logFEVrat_permVars.RData")

#Load the above saved data in the current workspace
load("nsub-p_REP_logFEVrat_permVars.RData")

#Gather the "gold-standard"/permutation-based FWERs, MWSLs, and Meffs
#Initialize empty vectors (all NA elements) for permutation-based gold-standard FWERs, MWSLs, and M-effs 
#for "identical" and multivariate-Normal transformed metabolites

PermFWER_Gold_iden.nsub = PermFWER_Gold_mN.nsub = 
  PermMWSL_Gold_iden.nsub = PermMWSL_Gold_mN.nsub = 
  PermMeff_Gold_iden.nsub = PermMeff_Gold_mN.nsub = rep(NA, REP)

#initiate null-vector to store the time taken for the gold-standard computations performed in the cluster computing system
#RER = R repetitions done in R parallel computations, nPerm = 10K (# of permutations)
timeT.nsub.10K = c()  
for(k in 1:REP){
  file_name = paste0("Part", k, ".RData")
  load(file_name)
  timeT.nsub.10K[k] = tocP[3]
  #for(j in 1:5){
  
  #print(c(k, j))
  print(k)
  
  #FWER
  PermFWER_Gold_iden.nsub[k] = rcopd_identity$t1err.percent/100
  PermFWER_Gold_mN.nsub[k] = rcopd_mn$t1err.percent/100
  
  #MWSL
  PermMWSL_Gold_iden.nsub[k] = rcopd_identity$res[1,1]
  PermMWSL_Gold_mN.nsub[k] = rcopd_mn$res[1,1]
  
  #M-Effs
  PermMeff_Gold_iden.nsub[k] = rcopd_identity$res[1,4]
  PermMeff_Gold_mN.nsub[k] = rcopd_mn$res[1,4]
  # }
}

mean(timeT.nsub.10K) 
#On an average, XXXXX secs taken for computing 10K permutation-based "identical" and "mN" versions

#Closed form expression: Meff using empirical Pearson correlation matrices for p features
#Initialize matrices for storing empirical correlation-based Meffs for different eigen-analysis-based Meff estimating methods
mat.rcopd_Meff_ecorr_sub.nsub = mat.rcopd_Meff_spcorr_sub.nsub = matrix(NA, REP, 7)
#Name the columns as per the methods
colnames(mat.rcopd_Meff_ecorr_sub.nsub) = colnames(mat.rcopd_Meff_spcorr_sub.nsub) =
  c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021")

source('Meff_sp.R') #connect the 'Meff_sp()' function in the current session 

for(k in 1:REP){
  print(k)
  filename = paste(c("logfeat_sub_res_nsub-p_", k, ".csv"), collapse = "")
  logfeat_sub_res = read.table(filename, sep = ",")
  feat = as.data.frame(logfeat_sub_res)
  
  #Pearson correlation
  rcopd_Meff_ecorr <- Meff(features=feat,
                           n.permutation=10^5,
                           method='ecorr',
                           big.mat=FALSE,
                           alpha=0.05)
  
  mat.rcopd_Meff_ecorr_sub.nsub[k,] <- c(
    rcopd_Meff_ecorr$Meff_Bonferroni,
    rcopd_Meff_ecorr$Meff_Sidak,
    rcopd_Meff_ecorr$Meff_Nyholt,
    rcopd_Meff_ecorr$Meff_Liji,
    rcopd_Meff_ecorr$Meff_Gao,
    rcopd_Meff_ecorr$Meff_Galwey,
    rcopd_Meff_ecorr$Meff_MWSL)
  
  #Spearman correlation
  rcopd_Meff_spcorr <- Meff_sp(features=feat,
                               n.permutation=10^5,
                               method='spcorr',
                               big.mat=FALSE,
                               alpha=0.05)
  
  mat.rcopd_Meff_spcorr_sub.nsub[k,] <- c(
    rcopd_Meff_spcorr$Meff_Bonferroni,
    rcopd_Meff_spcorr$Meff_Sidak,
    rcopd_Meff_spcorr$Meff_Nyholt,
    rcopd_Meff_spcorr$Meff_Liji,
    rcopd_Meff_spcorr$Meff_Gao,
    rcopd_Meff_spcorr$Meff_Galwey,
    rcopd_Meff_spcorr$Meff_MWSL)
}



#Distance correlation-based calculations
#install.packages("matrixcalc")
library(matrixcalc) #load the R-package "matrixcalc", for the "is.positive.definite()" function used below

liji_fun = function(x){
  val = (x>=1) + (x-floor(x))
  return(val)
}

gao_fun = function(x){
  numel = length(x)
  propel = rep(NA, numel)
  for(j in 1:numel){
    propel[j] = sum(x[1:j])/sum(x)
  }
  index = min(which(propel >= 0.995))
  return(index)
}

#Initialize NA-vectors, each of length REP, to store the M-effs obtained from distance correlation-based computations
Meff_dcorrNyholt.nsub = Meff_dcorrLiJi.nsub = Meff_dcorrGao.nsub = 
  Meff_dcorrGalwey.nsub = Meff_dcorrPeluso.nsub = rep(NA, REP)

#Check if the nearest-SPD matrices are indeed so, and assign eigen-values less than 10^{-12} to 0
#Distance correlation matrices are computed in MATLAB (see the MATLAb script "COPDGene_dcorr_logfeatures_p-nsub_logFEVrat_REP.m") 
checkPD = matrix(NA, REP, 4)

#Check if the sum of eigen-values exceed the total data variability
sum.eigvals.dcorr.1 = sum.eigvals.dcorr.2 = rep(NA, REP)

for(k in 1:REP){
  print(k)
  filename = paste(c("COPDGene_nearPD_dcorr_logfeatures_sub_res_p-nsub_logFEVrat_", k, ".txt"),
                   collapse = "")
  nearPD_dcorr_logfeatures = read.table(filename, header = F, sep=",")
  nearPD_dcorr_logfeatures = as.matrix(nearPD_dcorr_logfeatures)
  checkPD[k,1] = is.positive.definite(nearPD_dcorr_logfeatures) #check for PD
  eig.dcorr = eigen(nearPD_dcorr_logfeatures)
  eig.vals.dcorr = eig.dcorr$values
  
  sum.eigvals.dcorr.1[k] = ifelse(sum(eig.vals.dcorr) > sum(diag(nearPD_dcorr_logfeatures)), 1, 0)
  
  checkPD[k,2] = sum(eig.vals.dcorr < 0) #check total number of eigen-values less than 0
  checkPD[k,3] = sum(eig.vals.dcorr <= 10^(-12)) #check total number of eigen-values less than 10^{-12}
  if(!checkPD[k,1] & checkPD[k,2] > 0){
    checkPD[k,4] = min(eig.vals.dcorr[(eig.vals.dcorr < 0)]) #check the minimum of the negative eigen-values
  }else{checkPD[k,4] = 0}
  
  eig.vals.dcorr[(eig.vals.dcorr < 10^(-12))] = 0 #assign eigen-values smaller than 10^{-12} to 0
  sum.eigvals.dcorr.2[k] = ifelse(sum(eig.vals.dcorr) > sum(diag(nearPD_dcorr_logfeatures)), 1, 0)
  
  
  M = ncol(nearPD_dcorr_logfeatures)
  #Meff_Nyholt
  Meff_dcorrNyholt.nsub[k] = 1+(M-1)*(1-var(eig.vals.dcorr)/M)
  
  #Meff_JiLi
  Meff_dcorrLiJi.nsub[k] = sum(liji_fun(abs(eig.vals.dcorr)))
  
  #Meff_Gao
  Meff_dcorrGao.nsub[k] = gao_fun(eig.vals.dcorr)
  
  #Meff_Galwey
  Meff_dcorrGalwey.nsub[k] = sum(sqrt(eig.vals.dcorr))^2/sum(eig.vals.dcorr)
  
  #Meff_dcorr2021
  Meff_dcorrPeluso.nsub[k] = (sum(sqrt(eig.vals.dcorr))/log(eig.vals.dcorr[1]))^2/(sum(eig.vals.dcorr)/eig.vals.dcorr[1] + sqrt(eig.vals.dcorr[1]))
}

#Check sum of eigen-values
#raw eigen-analysis on near_PD matrices
summary(sum.eigvals.dcorr.1)
table(sum.eigvals.dcorr.1)
#after eigen-values < 10^{-12} assigned to 0 
summary(sum.eigvals.dcorr.2) #all 0's
table(sum.eigvals.dcorr.2) #all 0's

#check the distribution of the eigen-values that are minimum among the negative ones in each REP
summary(checkPD[,4]) #essentially zeroes, all good!


############## Use grouped information for p=761 metabolites into 8 Super Pathways #####################
#Closed form expression: Meff using empirical Pearson correlation matrices for p=761 features
########################################################################################################

#Empirical (Pearson) correlation-based Meff calculations ##############
#Initialize NA arrays of dimensions (REP, 7, 8) for both PrsCo and DisCo-based calculations
mat.rcopd_Meff_ecorr_sub_gp.nsub = mat.rcopd_Meff_spcorr_sub_gp.nsub = array(NA, dim = c(REP, 7, 8))
mat.rcopd_Meff_ecorr_sub_gpFin.nsub = mat.rcopd_Meff_spcorr_sub_gpFin.nsub = matrix(NA, REP, 7)

#name the array-dimensions - REPs, methods, super-pathway-based groups
dimnames(mat.rcopd_Meff_ecorr_sub_gp.nsub) = dimnames(mat.rcopd_Meff_spcorr_sub_gp.nsub) =
  list(as.character(1:REP), c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021"),
       names(sup.path.ind.v2))
colnames(mat.rcopd_Meff_ecorr_sub_gpFin.n100) = colnames(mat.rcopd_Meff_spcorr_sub_gpFin.n100) =
  c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021")

for(k in 1:REP){
  print(k)
  for(l in 1:8){
    filename = paste(c('logmet_class', l, '_sub_res[[', k, ']]'), collapse = "")
    feat = as.data.frame(eval(parse(text = filename)))
    rcopd_Meff_ecorr <- Meff(features=feat,
                             n.permutation=10^5,
                             method='ecorr',
                             big.mat=FALSE,
                             alpha=0.05)
    
    mat.rcopd_Meff_ecorr_sub_gp.nsub[k,,l] <- c(
      rcopd_Meff_ecorr$Meff_Bonferroni,
      rcopd_Meff_ecorr$Meff_Sidak,
      rcopd_Meff_ecorr$Meff_Nyholt,
      rcopd_Meff_ecorr$Meff_Liji,
      rcopd_Meff_ecorr$Meff_Gao,
      rcopd_Meff_ecorr$Meff_Galwey,
      rcopd_Meff_ecorr$Meff_MWSL)
    
    rcopd_Meff_spcorr <- Meff_sp(features=feat,
                                 n.permutation=10^5,
                                 method='spcorr',
                                 big.mat=FALSE,
                                 alpha=0.05)
    
    mat.rcopd_Meff_spcorr_sub_gp.nsub[k,,l] <- c(
      rcopd_Meff_spcorr$Meff_Bonferroni,
      rcopd_Meff_spcorr$Meff_Sidak,
      rcopd_Meff_spcorr$Meff_Nyholt,
      rcopd_Meff_spcorr$Meff_Liji,
      rcopd_Meff_spcorr$Meff_Gao,
      rcopd_Meff_spcorr$Meff_Galwey,
      rcopd_Meff_spcorr$Meff_MWSL)
  }
  
  for(methods in 1:7)
  {
    mat.rcopd_Meff_ecorr_sub_gpFin.nsub[k,methods] = sum(mat.rcopd_Meff_ecorr_sub_gp.nsub[k,methods,])
    mat.rcopd_Meff_spcorr_sub_gpFin.nsub[k,methods] = sum(mat.rcopd_Meff_spcorr_sub_gp.nsub[k,methods,])
  }
}

###################### Distance correlation based Meff calculations #######################

#Initiate NA matrices of dimentions REP-by-8 to store the Meffs obtained from distance correlation based calculations 
#corresponding to the 8 super-pathway-based groups
Meff_dcorrNyholt_gp.nsub = Meff_dcorrLiJi_gp.nsub = Meff_dcorrGao_gp.nsub = 
  Meff_dcorrGalwey_gp.nsub = Meff_dcorrPeluso_gp.nsub = matrix(NA, REP, 8)

#Check for PD and assign eigen-values less than 10^{-12} to 0
checkPD = array(NA, dim = c(REP, 4, 8))

#Check if the sum of eigen-values exceed M
sum.eigvals.dcorr.def.1 = sum.eigvals.dcorr.def.2 = matrix(0, REP, 8)
diff.eigvals.dcorr.def.1 = diff.eigvals.dcorr.def.2 = matrix(NA, REP, 8)


for(k in 1:REP){
  print(k)
  spathsize = c(365,99,181,24,25,32,10,25)
  for(l in 1:8){
    filename = paste(c("COPDGene_nearPD_dcorr_logfeatures_p",spathsize[l], "nsub_class", l,"_logFEVrat_", k, ".txt"),
                     collapse = "")
    nearPD_dcorr_logfeatures = read.table(filename, header = F, sep=",")
    nearPD_dcorr_logfeatures = as.matrix(nearPD_dcorr_logfeatures)
    checkPD[k,1,l] = is.positive.definite(nearPD_dcorr_logfeatures)
    eig.dcorr = eigen(nearPD_dcorr_logfeatures)
    eig.vals.dcorr = eig.dcorr$values
    
    if(sum(eig.vals.dcorr) > sum(diag(nearPD_dcorr_logfeatures))){
      sum.eigvals.dcorr.def.1[k,l] = 1
      diff.eigvals.dcorr.def.1[k, l] = abs(sum(eig.vals.dcorr) - sum(diag(nearPD_dcorr_logfeatures)))
    }
    
    checkPD[k,2,l] = sum(eig.vals.dcorr < 0)
    checkPD[k,3,l] = sum(eig.vals.dcorr <= 10^(-12))
    if(!checkPD[k,1,l] & checkPD[k,2,l] > 0){
      checkPD[k,4,l] = min(eig.vals.dcorr[(eig.vals.dcorr < 0)])
    }else{checkPD[k,4,l] = 0}
    eig.vals.dcorr[(eig.vals.dcorr < 10^(-12))] = 0
    
    if(sum(eig.vals.dcorr) > sum(diag(nearPD_dcorr_logfeatures))){
      sum.eigvals.dcorr.def.2[k,l] = 1
      diff.eigvals.dcorr.def.2[k, l] = abs(sum(eig.vals.dcorr) - sum(diag(nearPD_dcorr_logfeatures)))
    }
    
    M = ncol(nearPD_dcorr_logfeatures)
    #Meff_Nyholt
    Meff_dcorrNyholt_gp.nsub[k,l] = 1+(M-1)*(1-var(eig.vals.dcorr)/M)
    
    #Meff_JiLi
    Meff_dcorrLiJi_gp.nsub[k,l] = sum(liji_fun(abs(eig.vals.dcorr)))
    
    #Meff_Gao
    Meff_dcorrGao_gp.nsub[k,l] = gao_fun(eig.vals.dcorr)
    
    #Meff_Galwey
    Meff_dcorrGalwey_gp.nsub[k,l] = sum(sqrt(eig.vals.dcorr))^2/sum(eig.vals.dcorr)
    
    #Meff_dcorr2021
    Meff_dcorrPeluso_gp.nsub[k,l] = (sum(sqrt(eig.vals.dcorr))/log(eig.vals.dcorr[1]))^2/(sum(eig.vals.dcorr)/eig.vals.dcorr[1] + sqrt(eig.vals.dcorr[1]))
  }
}

#Check sum of eigen-values
#raw eigen-analysis on near_PD matrices
apply(sum.eigvals.dcorr.def.1, 2, table)
apply(diff.eigvals.dcorr.def.1, 2, summary)
#after eigen-values < 10^{-12} assigned to 0 
apply(sum.eigvals.dcorr.def.2, 2, table)
apply(diff.eigvals.dcorr.def.2, 2, summary)


#Check the distribution of minimum of the negative eigen-values
for(l in 1:8){
  print(summary(checkPD[,4,l])) #essentially zeroes
}



############## Use RANDOM grouped information for p=p metabolites into 8 Super Pathways #####################
#Closed form expression: Meff using empirical Pearson correlation, Spearman correlation, and distance correlation matrices 
#for p=p features
##########################################################################################################################

#Empirical (Pearson) correlation-based Meff calculations ##############
#Initiate arrays of dimensions (REP, 7, 8) to store Meffs corresponding to both Pearson and Spearman correlation matrices
mat.rcopd_Meff_ecorr_sub_rangp.nsub = mat.rcopd_Meff_spcorr_sub_rangp.nsub = array(NA, dim = c(REP, 7, 8))
mat.rcopd_Meff_ecorr_sub_rangpFin.nsub = mat.rcopd_Meff_spcorr_sub_rangpFin.nsub = matrix(NA, REP, 7)

#name the array-dimensions, as per REPs, methods, and super-pathway-based metabolite groups
dimnames(mat.rcopd_Meff_ecorr_sub_rangp.nsub) = dimnames(mat.rcopd_Meff_spcorr_sub_rangp.nsub) = 
  list(as.character(1:REP), c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021"),
       names(sup.path.ind.v2))
colnames(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub) = colnames(mat.rcopd_Meff_spcorr_sub_rangpFin.nsub) =  
  c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021")

#Obtain the RANDOM ALLOCATION indices (generated in MATLAB)
rand0Ind_nsub-p_Final = as.numeric(read.table("rand0Ind_nsub-p_Final.txt",
                                                sep = ',', header = F))
class(rand0Ind_nsub-p_Final)
p_sub = length(rand0Ind_nsub-p_Final);
p_sub_8class = c(0,365, 99, 181, 24, 25, 32, 10, 25)
p_sub_8class_cumSum = cumsum(p_sub_8class)

for(k in 1:REP){
  print(k)
  for(l in 1:8){
    filename = paste(c("logfeat_sub_res_nsub-p_", k, ".csv"), collapse = "")
    logfeat_sub_res = read.table(filename, sep = ",")
    feat = as.data.frame(logfeat_sub_res)[,rand0Ind_nsub-p_Final[(p_sub_8class_cumSum[l]+1):p_sub_8class_cumSum[l+1]]]
    rcopd_Meff_ecorr <- Meff(features=feat,
                             n.permutation=10^5,
                             method='ecorr',
                             big.mat=FALSE,
                             alpha=0.05)
    
    mat.rcopd_Meff_ecorr_sub_rangp.nsub[k,,l] <- c(
      rcopd_Meff_ecorr$Meff_Bonferroni,
      rcopd_Meff_ecorr$Meff_Sidak,
      rcopd_Meff_ecorr$Meff_Nyholt,
      rcopd_Meff_ecorr$Meff_Liji,
      rcopd_Meff_ecorr$Meff_Gao,
      rcopd_Meff_ecorr$Meff_Galwey,
      rcopd_Meff_ecorr$Meff_MWSL)
    
    rcopd_Meff_spcorr <- Meff_sp(features=feat,
                                 n.permutation=10^5,
                                 method='spcorr',
                                 big.mat=FALSE,
                                 alpha=0.05)
    
    mat.rcopd_Meff_spcorr_sub_rangp.nsub[k,,l] <- c(
      rcopd_Meff_spcorr$Meff_Bonferroni,
      rcopd_Meff_spcorr$Meff_Sidak,
      rcopd_Meff_spcorr$Meff_Nyholt,
      rcopd_Meff_spcorr$Meff_Liji,
      rcopd_Meff_spcorr$Meff_Gao,
      rcopd_Meff_spcorr$Meff_Galwey,
      rcopd_Meff_spcorr$Meff_MWSL)
  }
  
  for(methods in 1:7)
  {
    mat.rcopd_Meff_ecorr_sub_rangpFin.nsub[k,methods] = sum(mat.rcopd_Meff_ecorr_sub_rangp.nsub[k,methods,])
    mat.rcopd_Meff_spcorr_sub_rangpFin.nsub[k,methods] = sum(mat.rcopd_Meff_spcorr_sub_rangp.nsub[k,methods,])
  }
}

###################### Distance correlation based Meff calculations #######################

#Store the results corresponding to Distance Correlation
#Initiate matrices of dimensions REP-by-8
Meff_dcorrNyholt_rangp.nsub = Meff_dcorrLiJi_rangp.nsub = Meff_dcorrGao_rangp.nsub = 
  Meff_dcorrGalwey_rangp.nsub = Meff_dcorrPeluso_rangp.nsub = matrix(NA, REP, 8)

#Check for PD and assign eigen-values less than 10^{-12} to 0
checkPD = array(NA, dim = c(REP, 4, 8))

#Check if the sum of eigen-values exceed M
sum.eigvals.dcorr.ran.1 = sum.eigvals.dcorr.ran.2 = matrix(0, REP, 8)
diff.eigvals.dcorr.ran.1 = diff.eigvals.dcorr.ran.2 = matrix(NA, REP, 8)

for(k in 1:REP){
  print(k)
  spathsize = c(365,99,181,24,25,32,10,25)
  for(l in 1:8){
    filename = paste(c("COPDGene_nearPD_dcorr_logfeatures_p",spathsize[l], "nsub_randclass", l,"_logFEVrat_rand0_", k, ".txt"),
                     collapse = "")
    nearPD_dcorr_logfeatures = read.table(filename, header = F, sep=",")
    nearPD_dcorr_logfeatures = as.matrix(nearPD_dcorr_logfeatures)
    checkPD[k,1,l] = is.positive.definite(nearPD_dcorr_logfeatures)
    eig.dcorr = eigen(nearPD_dcorr_logfeatures)
    eig.vals.dcorr = eig.dcorr$values
    
    if(sum(eig.vals.dcorr) > sum(diag(nearPD_dcorr_logfeatures))){
      sum.eigvals.dcorr.ran.1[k,l] = 1
      diff.eigvals.dcorr.ran.1[k, l] = abs(sum(eig.vals.dcorr) - sum(diag(nearPD_dcorr_logfeatures)))
    }
    
    checkPD[k,2,l] = sum(eig.vals.dcorr < 0)
    checkPD[k,3,l] = sum(eig.vals.dcorr <= 10^(-12))
    if(!checkPD[k,1,l] &  checkPD[k,2,l] > 0){
      checkPD[k,4,l] = min(eig.vals.dcorr[(eig.vals.dcorr < 10^(-12))])
    }else{checkPD[k,4,l] = 0}
    
    eig.vals.dcorr[(eig.vals.dcorr < 10^(-12))] = 0
    
    if(sum(eig.vals.dcorr) > sum(diag(nearPD_dcorr_logfeatures))){
      sum.eigvals.dcorr.ran.2[k,l] = 1
      diff.eigvals.dcorr.ran.2[k, l] = abs(sum(eig.vals.dcorr) - sum(diag(nearPD_dcorr_logfeatures)))
    }
    
    M = ncol(nearPD_dcorr_logfeatures)
    #Meff_Nyholt
    Meff_dcorrNyholt_rangp.nsub[k,l] = 1+(M-1)*(1-var(eig.vals.dcorr)/M)
    
    #Meff_JiLi
    Meff_dcorrLiJi_rangp.nsub[k,l] = sum(liji_fun(abs(eig.vals.dcorr)))
    
    #Meff_Gao
    Meff_dcorrGao_rangp.nsub[k,l] = gao_fun(eig.vals.dcorr)
    
    #Meff_Galwey
    Meff_dcorrGalwey_rangp.nsub[k,l] = sum(sqrt(eig.vals.dcorr))^2/sum(eig.vals.dcorr)
    
    #Meff_dcorr2021
    Meff_dcorrPeluso_rangp.nsub[k,l] = (sum(sqrt(eig.vals.dcorr))/log(eig.vals.dcorr[1]))^2/(sum(eig.vals.dcorr)/eig.vals.dcorr[1] + sqrt(eig.vals.dcorr[1]))
  }
}

#Check sum of eigen-values
#raw eigen-analysis on near_PD matrices
apply(sum.eigvals.dcorr.ran.1, 2, table)
apply(diff.eigvals.dcorr.ran.1, 2, summary)
#after eigen-values < 10^{-12} assigned to 0 
apply(sum.eigvals.dcorr.ran.2, 2, table)
apply(diff.eigvals.dcorr.ran.2, 2, summary)


#Check distribution of minimum eigen-values less than 0
for(l in 1:8){
  print(summary(checkPD[,4,l]))
}




#Save the results in an RData file
save(PermFWER_Gold_iden.nsub, PermFWER_Gold_mN.nsub, PermMWSL_Gold_iden.nsub,
     PermMWSL_Gold_mN.nsub, PermMeff_Gold_iden.nsub, PermMeff_Gold_mN.nsub,
     mat.rcopd_Meff_ecorr_sub.nsub, mat.rcopd_Meff_spcorr_sub.nsub, 
     Meff_dcorrNyholt.nsub, Meff_dcorrLiJi.nsub,
     Meff_dcorrGao.nsub, Meff_dcorrGalwey.nsub, Meff_dcorrPeluso.nsub,
     mat.rcopd_Meff_ecorr_sub_gp.nsub, mat.rcopd_Meff_spcorr_sub_gp.nsub,
     Meff_dcorrNyholt_gp.nsub, Meff_dcorrLiJi_gp.nsub,
     Meff_dcorrGao_gp.nsub, Meff_dcorrGalwey_gp.nsub, Meff_dcorrPeluso_gp.nsub,
     mat.rcopd_Meff_ecorr_sub_rangp.nsub, mat.rcopd_Meff_spcorr_sub_rangp.nsub,
     Meff_dcorrNyholt_rangp.nsub, Meff_dcorrLiJi_rangp.nsub,
     Meff_dcorrGao_rangp.nsub, Meff_dcorrGalwey_rangp.nsub, Meff_dcorrPeluso_rangp.nsub,
     file = "COPDGene_nsub-p_allResults.RData")



##############################################################################################
#############################         Check effect of RANDOM assignments of metabolites in 8 groups ###################
rm(list=ls())
REP = R
randRep = RR #no. of repetitions for the random metabolite grouping
load("final_processed_dataframes_subjects_metabolites.RData")
load("COPDGene_nsub-p_allResults.RData")
load("nsub-p_REP_logFEVrat_permVars.RData")

library(matrixcalc)

#Empirical (Pearson) correlation-based Meff calculations 
mat.rcopd_Meff_ecorr_sub_rangp.All.nsub = array(NA, dim = c(REP, 7, 8, randRep))
mat.rcopd_Meff_ecorr_sub_rangpFin.All.nsub = array(NA, dim = c(REP, 7, randRep))

dimnames(mat.rcopd_Meff_ecorr_sub_rangp.All.nsub) = 
  list(as.character(1:REP), c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021"),
       names(sup.path.ind.v2), as.character(paste0("randRep",1:randRep)))
dimnames(mat.rcopd_Meff_ecorr_sub_rangpFin.All.nsub) = 
  list(as.character(1:REP), c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021"),
       as.character(paste0("randRep",1:randRep)))


#Distance correlation based Meff calculations
Meff_dcorrNyholt_rangp.All.nsub = Meff_dcorrLiJi_rangp.All.nsub = Meff_dcorrGao_rangp.All.nsub = 
  Meff_dcorrGalwey_rangp.All.nsub = Meff_dcorrPeluso_rangp.All.nsub = array(NA, dim = c(REP, 8, randRep))

liji_fun = function(x){
  val = (x>=1) + (x-floor(x))
  return(val)
}

gao_fun = function(x){
  numel = length(x)
  propel = rep(NA, numel)
  for(j in 1:numel){
    propel[j] = sum(x[1:j])/sum(x)
  }
  index = min(which(propel >= 0.995))
  return(index)
}

p_sub_8class = c(0,365, 99, 181, 24, 25, 32, 10, 25)
p_sub_8class_cumSum = cumsum(p_sub_8class)

checkPD = array(NA, dim = c(REP, 4, 8, randRep))

tic = proc.time()
for(rRep in 1:randRep){
  #Obtain the RANDOM ALLOCATION indices (generated in MATLAB)
  filename = paste0("nsub-p_REP_logFEVrat/rand", rRep, "Ind_nsub-p_Final.txt")
  randInd_nsub-p_Final = as.numeric(read.table(filename, sep = ',', header = F))

  for(k in 1:REP){
    print(paste0("randRep : ", rRep, " PearCorr : ", k))
    for(l in 1:8){
      feat = as.data.frame(logfeat_sub_res[[k]])[,randInd_nsub-p_Final[(p_sub_8class_cumSum[l]+1):p_sub_8class_cumSum[l+1]]]
      rcopd_Meff_ecorr <- Meff(features=feat,
                               n.permutation=10^5,
                               method='ecorr',
                               big.mat=FALSE,
                               alpha=0.05)
      
      mat.rcopd_Meff_ecorr_sub_rangp.All.nsub[k,,l,rRep] <- c(
        rcopd_Meff_ecorr$Meff_Bonferroni,
        rcopd_Meff_ecorr$Meff_Sidak,
        rcopd_Meff_ecorr$Meff_Nyholt,
        rcopd_Meff_ecorr$Meff_Liji,
        rcopd_Meff_ecorr$Meff_Gao,
        rcopd_Meff_ecorr$Meff_Galwey,
        rcopd_Meff_ecorr$Meff_MWSL)
    }
    
    for(methods in 1:7)
    {
      mat.rcopd_Meff_ecorr_sub_rangpFin.All.nsub[k,methods,rRep] = sum(mat.rcopd_Meff_ecorr_sub_rangp.All.nsub[k,methods,,rRep])
    }
  }
  
  for(k in 1:REP){
    print(paste0("randRep : ", rRep, " DistCorr : ", k))
    spathsize = c(365,99,181,24,25,32,10,25)
    for(l in 1:8){
      filename = paste(c("COPDGene_nearPD_dcorr_logfeatures_p",spathsize[l], "nsub_randclass", l,"_logFEVrat_rand", rRep, "_", k, ".txt"),
                       collapse = "")
      nearPD_dcorr_logfeatures = read.table(filename, header = F, sep=",")
      nearPD_dcorr_logfeatures = as.matrix(nearPD_dcorr_logfeatures)
      checkPD[k,1,l,rRep] = is.positive.definite(nearPD_dcorr_logfeatures)
      eig.dcorr = eigen(nearPD_dcorr_logfeatures)
      eig.vals.dcorr = eig.dcorr$values
      checkPD[k,2,l,rRep] = sum(eig.vals.dcorr < 0)
      checkPD[k,3,l,rRep] = sum(eig.vals.dcorr <= 10^(-12))
      if(!checkPD[k,1,l,rRep]){
        checkPD[k,4,l,rRep] = min(eig.vals.dcorr[(eig.vals.dcorr < 10^(-12))])
      }
      eig.vals.dcorr[(eig.vals.dcorr < 10^(-12))] = 0
      
      M = ncol(nearPD_dcorr_logfeatures)
      #Meff_Nyholt
      Meff_dcorrNyholt_rangp.All.nsub[k,l,rRep] = 1+(M-1)*(1-var(eig.vals.dcorr)/M)
      
      #Meff_JiLi
      Meff_dcorrLiJi_rangp.All.nsub[k,l,rRep] = sum(liji_fun(abs(eig.vals.dcorr)))
      
      #Meff_Gao
      Meff_dcorrGao_rangp.All.nsub[k,l,rRep] = gao_fun(eig.vals.dcorr)
      
      #Meff_Galwey
      Meff_dcorrGalwey_rangp.All.nsub[k,l,rRep] = sum(sqrt(eig.vals.dcorr))^2/sum(eig.vals.dcorr)
      
      #Meff_dcorr2021
      Meff_dcorrPeluso_rangp.All.nsub[k,l,rRep] = (sum(sqrt(eig.vals.dcorr))/log(eig.vals.dcorr[1]))^2/(sum(eig.vals.dcorr)/eig.vals.dcorr[1] + sqrt(eig.vals.dcorr[1]))
    }
  }
}
toc = proc.time() - tic
toc

#save the key data objects
save(mat.rcopd_Meff_ecorr_sub_rangp.All.nsub, mat.rcopd_Meff_ecorr_sub_rangpFin.All.nsub,
     Meff_dcorrNyholt_rangp.All.nsub, Meff_dcorrLiJi_rangp.All.nsub, Meff_dcorrGao_rangp.All.nsub,
     Meff_dcorrGalwey_rangp.All.nsub, Meff_dcorrPeluso_rangp.All.nsub,
     file = "COPDGene_nsub-p_randRep.RData")

#Summarize the random replications results
rm(list=ls())
load("COPDGene_nsub-p_randRep.RData")

REP = R
randRep = RR #no. of repetitions for random assignment


mat.rcopd_Meff_dcorr_sub_rangpFin.Nyholt.nsub = mat.rcopd_Meff_dcorr_sub_rangpFin.LiJi.nsub =
  mat.rcopd_Meff_dcorr_sub_rangpFin.Gao.nsub = mat.rcopd_Meff_dcorr_sub_rangpFin.Galwey.nsub = 
  mat.rcopd_Meff_dcorr_sub_rangpFin.Peluso.nsub = array(NA, dim = c(REP, randRep))


for(rep in 1:REP){
  for(rand in 1:randRep){
    mat.rcopd_Meff_dcorr_sub_rangpFin.Nyholt.nsub[rep, rand] = sum(Meff_dcorrNyholt_rangp.All.nsub[rep, , rand])
    mat.rcopd_Meff_dcorr_sub_rangpFin.LiJi.nsub[rep, rand] = sum(Meff_dcorrLiJi_rangp.All.nsub[rep, , rand])
    mat.rcopd_Meff_dcorr_sub_rangpFin.Gao.nsub[rep, rand] = sum(Meff_dcorrGao_rangp.All.nsub[rep, , rand])
    mat.rcopd_Meff_dcorr_sub_rangpFin.Galwey.nsub[rep, rand] = sum(Meff_dcorrGalwey_rangp.All.nsub[rep, , rand])
    mat.rcopd_Meff_dcorr_sub_rangpFin.Peluso.nsub[rep, rand] = sum(Meff_dcorrPeluso_rangp.All.nsub[rep, , rand])
  }
}

class(mat.rcopd_Meff_ecorr_sub_rangpFin.All.nsub)
dim(mat.rcopd_Meff_ecorr_sub_rangpFin.All.nsub) #REP: methods: randRep = R x 7 x RR

#save density plots for each randRep, each method
denREP.methods.randRep = array(list(), dim = c(7, randRep)) 
#save 5-number-summary for each randRep, each method
methodNames = c("Bonferroni", "Sidak", "Nyholt", "LiJi", "Gao", "Galwey", "Peluso")

fiveNum.methods.randRep = array(NA, dim = c(6, 7, RR), 
                                dimnames = list(c("Min","Q1", "Median", "Mean", "Q3", "Max"), 
                                                methodNames,
                                                paste(1:RR))
)


for(method in 3:7){
  
  filename1 = paste(c("randRep_PrsCo_", methodNames[method], "_bp_randGrp_logFEVrat_nsub-p_REP_mN.png"),
                    collapse = "")
  png(filename1, width = 20, height = 12, units = "in", res = 300)
  boxplot(mat.rcopd_Meff_ecorr_sub_rangpFin.All.nsub[, method, ], 
          main = paste(c("Boxplots of Estimated PrsCo M-effs | Repeated Random Metabolites Grouping | Method: ", methodNames[method], "\n n = nsub, Outcome = log-FEV-ratio"), collapse = ""),
          xlab = "Random Grouping Repetition #",
          ylab = "Distribution of Estimated M-eff's",
          col = 1:randRep,
          cex.lab = 1.25)
  dev.off()
  
  for(rand in 1:randRep){
    
    print(paste(c("Method: ", method, " | Rand:", rand), collapse = ""))
    fiveNum.methods.randRep[, method, rand] = summary(mat.rcopd_Meff_ecorr_sub_rangpFin.All.nsub[, method, rand])
    
  }
}

#DisCo: Nyholt
filename = paste(c("randRep_DisCo_Nyholt_bp_randGrp_logFEVrat_nsub-p_REP_mN.png"),
                 collapse = "")
png(filename, width = 20, height = 12, units = "in", res = 300)
boxplot(mat.rcopd_Meff_dcorr_sub_rangpFin.Nyholt.nsub, 
        main = "Boxplots of Estimated DisCo M-effs | Repeated Random Metabolites Grouping | Method: Nyholt \n n = nsub, Outcome: log-FEV-ratio",
        xlab = "Random Grouping Repetition #",
        ylab = "Distribution of Estimated M-eff's",
        col = 1:randRep,
        cex.lab = 1.25)
dev.off()

#DisCo: LiJi
filename = paste(c("randRep_DisCo_LiJi_bp_randGrp_logFEVrat_nsub-p_REP_mN.png"),
                 collapse = "")
png(filename, width = 20, height = 12, units = "in", res = 300)
boxplot(mat.rcopd_Meff_dcorr_sub_rangpFin.LiJi.n100, 
        main = "Boxplots of Estimated DisCo M-effs | Repeated Random Metabolites Grouping | Method: LiJi \n n = nsub, Outcome: log-FEV-ratio",
        xlab = "Random Grouping Repetition #",
        ylab = "Distribution of Estimated M-eff's",
        col = 1:randRep,
        cex.lab = 1.25)
dev.off()

#DisCo: Gao
filename = paste(c("randRep_DisCo_Gao_bp_randGrp_logFEVrat_nsub-p_REP_mN.png"),
                 collapse = "")
png(filename, width = 20, height = 12, units = "in", res = 300)
boxplot(mat.rcopd_Meff_dcorr_sub_rangpFin.Gao.n100, 
        main = "Boxplots of Estimated DisCo M-effs | Repeated Random Metabolites Grouping | Method: Gao \n n = nsub, Outcome: log-FEV-ratio",
        xlab = "Random Grouping Repetition #",
        ylab = "Distribution of Estimated M-eff's",
        col = 1:randRep,
        cex.lab = 1.25)
dev.off()

#DisCo: Galwey
filename = paste(c("randRep_DisCo_Galwey_bp_randGrp_logFEVrat_nsub-p_REP_mN.png"),
                 collapse = "")
png(filename, width = 20, height = 12, units = "in", res = 300)
boxplot(mat.rcopd_Meff_dcorr_sub_rangpFin.Galwey.n100, 
        main = "Boxplots of Estimated DisCo M-effs | Repeated Random Metabolites Grouping | Method: Galwey \n n = nsub, Outcome: log-FEV-ratio",
        xlab = "Random Grouping Repetition #",
        ylab = "Distribution of Estimated M-eff's",
        col = 1:randRep,
        cex.lab = 1.25)
dev.off()

#DisCo: Peluso
filename = paste(c("randRep_DisCo_Peluso_bp_randGrp_logFEVrat_nsub-p_REP_mN.png"),
                 collapse = "")
png(filename, width = 20, height = 12, units = "in", res = 300)
boxplot(mat.rcopd_Meff_dcorr_sub_rangpFin.Peluso.n100, 
        main = "Boxplots of Estimated DisCo M-effs | Repeated Random Metabolites Grouping | Method: Peluso \n n = nsub, Outcome: log-FEV-ratio",
        xlab = "Random Grouping Repetition #",
        ylab = "Distribution of Estimated M-eff's",
        col = 1:randRep,
        cex.lab = 1.25)
dev.off()




############################################################################
############# Plot the boxplots ############################################

#Metabolites Grouping: Ungrouped

#M-effs 
Meffs_comp_Gold.nsub = stack(data.frame(Perm_id = PermMeff_Gold_iden.nsub,
                                        Perm_mN = PermMeff_Gold_mN.nsub,
                                        Bon = mat.rcopd_Meff_ecorr_sub.nsub[,1],
                                        Sidak = mat.rcopd_Meff_ecorr_sub.nsub[,2],
                                        Nyholt = mat.rcopd_Meff_ecorr_sub.nsub[,3],
                                        dcNyholt = Meff_dcorrNyholt.nsub,
                                        LiJi = mat.rcopd_Meff_ecorr_sub.nsub[,4],
                                        dcLiJi = Meff_dcorrLiJi.nsub,
                                        Gao = mat.rcopd_Meff_ecorr_sub.nsub[,5],
                                        dcGao = Meff_dcorrGao.nsub,
                                        Galwey = mat.rcopd_Meff_ecorr_sub.nsub[,6],
                                        dcGalwey = Meff_dcorrGalwey.nsub,
                                        Peluso = mat.rcopd_Meff_ecorr_sub.nsub[,7],
                                        dcPeluso = Meff_dcorrPeluso.nsub
))

maxAll = max(Meffs_comp_Gold.nsub[,1])
minAll = min(Meffs_comp_Gold.nsub[,1])

png("allMeffs_logFEVrat_nsub-p_REP_nPerm10Kalpha05.png",
    width = 12, height = 8, units = "in", res = 300)
library(ggplot2)
bp.Meffs_comp_Gold = ggplot(Meffs_comp_Gold.nsub, aes(x = ind, y = values)) + 
  #geom_violin(width=1.4) + 
  geom_boxplot(color = c(2,2,3,4,5,5,6,6,7,7,8,8,9,9), alpha=0.5, outlier.color = 10, notch = T) +
  geom_hline(yintercept = 761, colour = 'red') +
  scale_y_continuous(breaks = c(seq(30, maxAll, 25)), sec.axis = sec_axis(trans=~.*1, breaks = c(seq(30, maxAll, 25))))
#geom_abline(intercept = 26/45) 
bp.Meffs_comp_Gold = bp.Meffs_comp_Gold + labs(y="Estimated M-effs (Initial p = p)", x = "Methods: PearCorr & DistCorr for eigen analysis-based methods") +
  labs(title = "Estimated Effective Number of Tests | Response: log-FEV-ratio\nInitial p = p | REP subsets of sample size n = nsub\nnPermutations = 10K | FWER (alpha) = 0.05")
bp.Meffs_comp_Gold +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face = "bold")) +
  theme(axis.title = element_text(face = "bold"))
dev.off()

#Gold standard: Identity
MWSL_comp_Gold_iden.nsub = stack(data.frame(Perm = PermFWER_Gold_iden.nsub,
                                            Bon = PermMWSL_Gold_iden.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,1],
                                            Sidak = PermMWSL_Gold_iden.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,2],
                                            Nyholt = PermMWSL_Gold_iden.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,3],
                                            dcNyholt = PermMWSL_Gold_iden.nsub*Meff_dcorrNyholt.nsub,
                                            LiJi = PermMWSL_Gold_iden.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,4],
                                            dcLiJi = PermMWSL_Gold_iden.nsub*Meff_dcorrLiJi.nsub,
                                            Gao = PermMWSL_Gold_iden.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,5],
                                            dcGao = PermMWSL_Gold_iden.nsub*Meff_dcorrGao.nsub,
                                            Galwey = PermMWSL_Gold_iden.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,6],
                                            dcGalwey = PermMWSL_Gold_iden.nsub*Meff_dcorrGalwey.nsub,
                                            Peluso = PermMWSL_Gold_iden.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,7],
                                            dcPeluso = PermMWSL_Gold_iden.nsub*Meff_dcorrPeluso.nsub
))

maxAll = round(max(MWSL_comp_Gold_iden.nsub[,1]),3)
minAll = round(min(MWSL_comp_Gold_iden.nsub[,1]),3)

png("allFWERs_logFEVrat_nsub-p_REP_nPerm10Kalpha05_iden.png",
    width = 12, height = 8, units = "in", res = 300)
library(ggplot2)
bp.MWSL_comp_Gold_iden = ggplot(MWSL_comp_Gold_iden.nsub, aes(x = ind, y = values)) + 
  #geom_violin(width=1.4) + 
  geom_boxplot(color = c(2,3,4,5,5,6,6,7,7,8,8,9,9), alpha=0.5, outlier.color = 10, notch = T) +
  #coord_cartesian(ylim=c(0.45,maxAll)) 
  scale_y_continuous(breaks = c(seq(0, maxAll, 0.0025)), sec.axis = sec_axis(trans=~.*1, breaks = c(seq(0, maxAll, 0.0025)))) +
  geom_hline(yintercept = 0.05, colour = "red") 
bp.MWSL_comp_Gold_iden = bp.MWSL_comp_Gold_iden + labs(y="Estimated FWERs (alpha = 0.05)", x = "Methods: PearCorr & DistCorr for eigen analysis-based methods") +
  labs(title = "Estimated Family-Wise Error RatesResponse: log-FEV-ratio\nInitial p = p, REP subsets of sample size n = nsub\nnPermutations = 10K, alpha = 0.05, met-transform = Identity")
bp.MWSL_comp_Gold_iden +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face = "bold")) +
  theme(axis.title = element_text(face = "bold"))
dev.off()

#Gold standard: mN
MWSL_comp_Gold_mN.nsub = stack(data.frame(
  Perm = PermFWER_Gold_mN.nsub,
  Bon = PermMWSL_Gold_mN.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,1],
  Sidak = PermMWSL_Gold_mN.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,2],
  Nyholt = PermMWSL_Gold_mN.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,3],
  dcNyholt = PermMWSL_Gold_mN.nsub*Meff_dcorrNyholt.nsub,
  LiJi = PermMWSL_Gold_mN.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,4],
  dcLiJi = PermMWSL_Gold_mN.nsub*Meff_dcorrLiJi.nsub,
  Gao = PermMWSL_Gold_mN.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,5],
  dcGao = PermMWSL_Gold_mN.nsub*Meff_dcorrGao.nsub,
  Galwey = PermMWSL_Gold_mN.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,6],
  dcGalwey = PermMWSL_Gold_mN.nsub*Meff_dcorrGalwey.nsub,
  Peluso = PermMWSL_Gold_mN.nsub*mat.rcopd_Meff_ecorr_sub.nsub[,7],
  dcPeluso = PermMWSL_Gold_mN.nsub*Meff_dcorrPeluso.nsub
))

maxAll = max(MWSL_comp_Gold_mN.nsub[,1])
minAll = min(MWSL_comp_Gold_mN.nsub[,1])

png("allFWERs_logFEVrat_nsub-p_REP_nPerm10Kalpha05_mN.png",
    width = 12, height = 8, units = "in", res = 300)
library(ggplot2)
bp.MWSL_comp_Gold_mN = ggplot(MWSL_comp_Gold_mN.nsub, aes(x = ind, y = values)) + 
  #geom_violin(width=1.4) + 
  geom_boxplot(color = c(2,3,4,5,5,6,6,7,7,8,8,9,9), alpha=0.5, outlier.color = 10, notch = T) +
  geom_hline(yintercept = 0.05, colour = 'red') +
  #coord_cartesian(ylim=c(minAll,maxAll)) 
  scale_y_continuous(breaks = c(seq(0, maxAll, 0.0025)), sec.axis = sec_axis(trans = ~.*1, breaks = c(seq(0, maxAll, 0.0025))))
#geom_abline(intercept = 26/45) 
bp.MWSL_comp_Gold_mN = bp.MWSL_comp_Gold_mN + labs(y="Estimated FWERs (alpha = 0.05)", x = "Methods: PearCorr & DistCorr for eigen analysis-based methods") +
  labs(title = "Estimated Family-Wise Error Rates | Response: log-FEV-ratio\nInitial p = p, REP subsets of sample size n = nsub\nnPermutations = 10K, alpha = 0.05, met-transform = MV-Normal")
bp.MWSL_comp_Gold_mN +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face = "bold")) +
  theme(axis.title = element_text(face = "bold"))
dev.off()

#########################################################################
# Similarly, plots can be made for Metabolites' Grouping = Defined and Metabolites' Grouping = Random
#########################################################################


