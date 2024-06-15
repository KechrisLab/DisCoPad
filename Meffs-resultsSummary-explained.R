#Compute medians of Meffs/FWERs for all the methods
#COPDGene, n = 100

rm(list=ls()) #clear all objects from session
load("COPDGene_nsub-p_allResults.RData")
source('rmse.R') #link the rmse.R function in the session

#Compute median and median absolute deviation (MAD) 

round(median(PermMeff_Gold_mN.nsub)) 
round(mad(PermMeff_Gold_mN.nsub)) 


#Ungrouped metabolites - compute medians and RMSE's of estimated M-eff's
#Pearson correlation
dim(mat.rcopd_Meff_ecorr_sub.nsub) 
colnames(mat.rcopd_Meff_ecorr_sub.nsub)

med.PrsCo.cont.unGrp.nsub = round(apply(mat.rcopd_Meff_ecorr_sub.nsub, 2, median))
rmse.PrsCo.cont.unGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_ecorr_sub.nsub))
names(rmse.PrsCo.cont.unGrp.nsub) = colnames(mat.rcopd_Meff_ecorr_sub.nsub)
rmse.PrsCo.cont.unGrp.nsub = round(apply(mat.rcopd_Meff_ecorr_sub.nsub, 2, rmse, true = PermMeff_Gold_mN.nsub))
rmse.PrsCo.cont.unGrp.nsub

cbind(med.PrsCo.cont.unGrp.nsub, rmse.PrsCo.cont.unGrp.nsub)

#Spearman correlation
dim(mat.rcopd_Meff_spcorr_sub.nsub) 
colnames(mat.rcopd_Meff_spcorr_sub.nsub)

med.SprCo.cont.unGrp.nsub = round(apply(mat.rcopd_Meff_spcorr_sub.nsub, 2, median))
rmse.SprCo.cont.unGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_spcorr_sub.nsub))
names(rmse.SprCo.cont.unGrp.nsub)= colnames(mat.rcopd_Meff_spcorr_sub.nsub)
rmse.SprCo.cont.unGrp.nsub = round(apply(mat.rcopd_Meff_spcorr_sub.nsub, 2, rmse, true = PermMeff_Gold_mN.nsub))
rmse.SprCo.cont.unGrp.nsub

cbind(med.SprCo.cont.unGrp.nsub, rmse.SprCo.cont.unGrp.nsub)

#DisCo
mat.rcopd_Meff_dcorr_sub.nsub = as.matrix(cbind(Meff_dcorrNyholt.nsub, Meff_dcorrLiJi.nsub, Meff_dcorrGao.nsub,
                                                Meff_dcorrGalwey.nsub, Meff_dcorrPeluso.nsub))
colnames(mat.rcopd_Meff_dcorr_sub.nsub) = colnames(mat.rcopd_Meff_ecorr_sub.nsub)[-c(1:2)]
colnames(mat.rcopd_Meff_dcorr_sub.nsub)

med.DisCo.cont.unGrp.nsub = round(apply(mat.rcopd_Meff_dcorr_sub.nsub, 2, median))
rmse.DisCo.cont.unGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_dcorr_sub.nsub))
names(rmse.DisCo.cont.unGrp.nsub) = colnames(mat.rcopd_Meff_dcorr_sub.nsub)
rmse.DisCo.cont.unGrp.nsub = round(apply(mat.rcopd_Meff_dcorr_sub.nsub, 2, rmse, true = PermMeff_Gold_mN.nsub))
rmse.DisCo.cont.unGrp.nsub

cbind(med.DisCo.cont.unGrp.nsub, rmse.DisCo.cont.unGrp.nsub)


#Grouped (defined) metabolites
#Pearson correlation
REP = R

dim(mat.rcopd_Meff_ecorr_sub_gp.nsub) #REP / 7 / 8

mat.rcopd_Meff_ecorr_sub_gpFin.nsub = matrix(NA, REP, 7)
colnames(mat.rcopd_Meff_ecorr_sub_gpFin.nsub) = c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021")

for(k in 1:REP){
  for(methods in 1:7){
    mat.rcopd_Meff_ecorr_sub_gpFin.nsub[k,methods] = sum(mat.rcopd_Meff_ecorr_sub_gp.nsub[k,methods,])
  }
}

colnames(mat.rcopd_Meff_ecorr_sub_gpFin.nsub)
round(apply(mat.rcopd_Meff_ecorr_sub_gpFin.nsub, 2, median))

med.PrsCo.cont.dfGrp.nsub = round(apply(mat.rcopd_Meff_ecorr_sub_gpFin.nsub, 2, median))
rmse.PrsCo.cont.dfGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_ecorr_sub_gpFin.nsub))
names(rmse.PrsCo.cont.dfGrp.nsub) = colnames(mat.rcopd_Meff_ecorr_sub_gpFin.nsub)
rmse.PrsCo.cont.dfGrp.nsub = round(apply(mat.rcopd_Meff_ecorr_sub_gpFin.nsub, 2, rmse,
                                         true = PermMeff_Gold_mN.nsub))
rmse.PrsCo.cont.dfGrp.nsub

cbind(med.PrsCo.cont.dfGrp.nsub, rmse.PrsCo.cont.dfGrp.nsub)

#Spearman Correlation
dim(mat.rcopd_Meff_spcorr_sub_gp.nsub) #R / 7 / 8

mat.rcopd_Meff_spcorr_sub_gpFin.nsub = matrix(NA, REP, 7)
colnames(mat.rcopd_Meff_spcorr_sub_gpFin.nsub) = c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021")

for(k in 1:REP){
  for(methods in 1:7){
    mat.rcopd_Meff_spcorr_sub_gpFin.nsub[k,methods] = sum(mat.rcopd_Meff_spcorr_sub_gp.nsub[k,methods,])
  }
}

colnames(mat.rcopd_Meff_spcorr_sub_gpFin.nsub)
round(apply(mat.rcopd_Meff_spcorr_sub_gpFin.nsub, 2, median))

med.SprCo.cont.dfGrp.nsub = round(apply(mat.rcopd_Meff_spcorr_sub_gpFin.nsub, 2, median))
rmse.SprCo.cont.dfGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_spcorr_sub_gpFin.nsub))
names(rmse.SprCo.cont.dfGrp.nsub) = colnames(mat.rcopd_Meff_spcorr_sub_gpFin.nsub)
rmse.SprCo.cont.dfGrp.nsub = round(apply(mat.rcopd_Meff_spcorr_sub_gpFin.nsub, 2, rmse,
                                         true = PermMeff_Gold_mN.nsub))
rmse.SprCo.cont.dfGrp.nsub

cbind(med.SprCo.cont.dfGrp.nsub, rmse.SprCo.cont.dfGrp.nsub)

#DisCo
mat.rcopd_Meff_dcorr_sub.dfGrp.nsub = as.matrix(cbind(rowSums(Meff_dcorrNyholt_gp.nsub), rowSums(Meff_dcorrLiJi_gp.nsub), 
                                                      rowSums(Meff_dcorrGao_gp.nsub), rowSums(Meff_dcorrGalwey_gp.nsub), rowSums(Meff_dcorrPeluso_gp.nsub)))
colnames(mat.rcopd_Meff_dcorr_sub.dfGrp.nsub) = colnames(mat.rcopd_Meff_ecorr_sub_gpFin.nsub)[-c(1:2)]
colnames(mat.rcopd_Meff_dcorr_sub.dfGrp.nsub)

med.DisCo.cont.dfGrp.nsub = round(apply(mat.rcopd_Meff_dcorr_sub.dfGrp.nsub, 2, median))
rmse.DisCo.cont.dfGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_dcorr_sub.dfGrp.nsub))
names(rmse.DisCo.cont.dfGrp.nsub) = colnames(mat.rcopd_Meff_dcorr_sub.dfGrp.nsub)
rmse.DisCo.cont.dfGrp.nsub = round(apply(mat.rcopd_Meff_dcorr_sub.dfGrp.nsub, 2, rmse, true = PermMeff_Gold_mN.nsub))
rmse.DisCo.cont.dfGrp.nsub

cbind(med.DisCo.cont.dfGrp.nsub, rmse.DisCo.cont.dfGrp.nsub)


#Grouped (random) metabolites
#Pearson correlation
REP = R

dim(mat.rcopd_Meff_ecorr_sub_rangp.nsub) #REP / 7 / 8

mat.rcopd_Meff_ecorr_sub_rangpFin.nsub = matrix(NA, REP, 7)
colnames(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub) = c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021")

for(k in 1:REP){
  for(methods in 1:7){
    mat.rcopd_Meff_ecorr_sub_rangpFin.nsub[k,methods] = sum(mat.rcopd_Meff_ecorr_sub_rangp.nsub[k,methods,])
  }
}

colnames(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub)
round(apply(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub, 2, median))

med.PrsCo.cont.rnGrp.nsub = round(apply(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub, 2, median))
rmse.PrsCo.cont.rnGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub))
names(rmse.PrsCo.cont.rnGrp.nsub) = colnames(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub)
rmse.PrsCo.cont.rnGrp.nsub = round(apply(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub, 2, rmse,
                                         true = PermMeff_Gold_mN.nsub))
rmse.PrsCo.cont.rnGrp.nsub

cbind(med.PrsCo.cont.rnGrp.nsub, rmse.PrsCo.cont.rnGrp.nsub)

#Spearman Correlation
mat.rcopd_Meff_spcorr_sub_rangpFin.nsub = matrix(NA, REP, 7)
colnames(mat.rcopd_Meff_spcorr_sub_rangpFin.nsub) = c("Meff_Bon", "Meff_Sidak", "Meff_Nyholt", "Meff_LiJi", "Meff_Gao", "Meff_Galwey", "Meff_2021")

for(k in 1:REP){
  for(methods in 1:7){
    mat.rcopd_Meff_spcorr_sub_rangpFin.nsub[k,methods] = sum(mat.rcopd_Meff_spcorr_sub_rangp.nsub[k,methods,])
  }
}

colnames(mat.rcopd_Meff_spcorr_sub_rangpFin.nsub)
round(apply(mat.rcopd_Meff_spcorr_sub_rangpFin.nsub, 2, median))

med.SprCo.cont.rnGrp.nsub = round(apply(mat.rcopd_Meff_spcorr_sub_rangpFin.nsub, 2, median))
rmse.SprCo.cont.rnGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_spcorr_sub_rangpFin.nsub))
names(rmse.SprCo.cont.rnGrp.nsub) = colnames(mat.rcopd_Meff_spcorr_sub_rangpFin.nsub)
rmse.SprCo.cont.rnGrp.nsub = round(apply(mat.rcopd_Meff_spcorr_sub_rangpFin.nsub, 2, rmse,
                                         true = PermMeff_Gold_mN.nsub))
rmse.SprCo.cont.rnGrp.nsub

cbind(med.SprCo.cont.rnGrp.nsub, rmse.SprCo.cont.rnGrp.nsub)


#DisCo
mat.rcopd_Meff_dcorr_sub.rnGrp.nsub = as.matrix(cbind(rowSums(Meff_dcorrNyholt_rangp.nsub), rowSums(Meff_dcorrLiJi_rangp.nsub), 
                                                      rowSums(Meff_dcorrGao_rangp.nsub), rowSums(Meff_dcorrGalwey_rangp.nsub), 
                                                      rowSums(Meff_dcorrPeluso_rangp.nsub)))
colnames(mat.rcopd_Meff_dcorr_sub.rnGrp.nsub) = colnames(mat.rcopd_Meff_ecorr_sub_rangpFin.nsub)[-c(1:2)]
colnames(mat.rcopd_Meff_dcorr_sub.rnGrp.nsub)

med.DisCo.cont.rnGrp.nsub = round(apply(mat.rcopd_Meff_dcorr_sub.rnGrp.nsub, 2, median))
rmse.DisCo.cont.rnGrp.nsub = rep(NA, ncol(mat.rcopd_Meff_dcorr_sub.rnGrp.nsub))
names(rmse.DisCo.cont.rnGrp.nsub) = colnames(mat.rcopd_Meff_dcorr_sub.rnGrp.nsub)
rmse.DisCo.cont.rnGrp.nsub = round(apply(mat.rcopd_Meff_dcorr_sub.rnGrp.nsub, 2, rmse, true = PermMeff_Gold_mN.nsub))
rmse.DisCo.cont.rnGrp.nsub

cbind(med.DisCo.cont.rnGrp.nsub, rmse.DisCo.cont.rnGrp.nsub)







