#COPD Gene 
#nsub = n, p=761, REP = R

#start with a clean session
rm(list=ls())

#call source function "FWERperm2.R" to compute perm-based metabolome-wide significance level (MWSL) and corresponding M-eff
source(file = 'FWERperm2.R')

#load the data-frames of covariates, features (metabolites), and outcomes
load("nsub-p_REP_logFEVrat_permVars.RData")

#load the essential libraries for parallel computing
library(doParallel)
library(doRNG)
library(corpcor)
library(MASS)
library(future.apply)

#part number
#If REP = 100 and no. of parts = 20, then in each part, no of repetitions = 100/20 = 5
part = 1

#5 replicates in each part
rep_ind = ((part-1)*5 + 1) : (part*5)

plan(multisession) ### Tell R to use multiple cores
ticP = proc.time()
output <- future_lapply(rep_ind, function(li){
  print(li)
  out = unlist(outcome_logFEVrat_sub[[li]])
  feat = as.data.frame(logfeat_sub[[li]])
  conf = as.data.frame(conf_sub[[li]])
  rcopd_identity = FWERperm2(outcome=log(out/(1-out)),
                             features=feat,
                             confounders=conf,
                             n.permutation=10000,
                             method="identity",
                             verbose=F,
                             seed = li)
  rcopd_mn = FWERperm2(outcome=log(out/(1-out)),
                       features=feat,
                       confounders=conf,
                       n.permutation=10000,
                       method="mN",
                       verbose=F,
                       seed = li)
  return(
    list(
      rcopd_identity = rcopd_identity,
      rcopd_mn = rcopd_mn
    )
  )
}, future.seed = part)
tocP = proc.time() - ticP

#name the output object elements as per the REP number  
names(output) = paste("REP", rep_ind, sep = "")
plan(sequential) ### Go back to single core

#save the results objects
filename=paste0("COPDGene_Meff_nsub-p_logFEVrat_nperm10K_part", part, ".RData")
save(output, tocP, rep_ind, file = filename)



