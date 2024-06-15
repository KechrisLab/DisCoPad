#Prepare the COPDGene Dataset for final analysis

#start with a clean R environment
rm(list=ls())

#load the imputed and normalized dataset
load("D:/OneDrive - The University of Colorado Denver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/OneDrive_2021-05-19/MetaboGuru_Data/COPDGene_Metabolon/COPDGene_Metabolon_Metabolite_MetaData12Sep2018.Rda")

#Collect the full metabolite data info; check column names and rownames
dat.met = metabo
dim(dat.met) #1392 x 13
colnames(dat.met)
head(rownames(dat.met))

#tabulate the metabolites based on super-pathways, including partially and non-annotated ones
table(dat.met$`SUPER PATHWAY`, useNA = "ifany")
#Check the number of super-pathways
length(table(dat.met$`SUPER PATHWAY`, useNA = "ifany")) #10
#check the range of metabolite numbers falling under different super-pathways
range(table(dat.met$`SUPER PATHWAY`, useNA = "ifany")) #(10, 431)

#Check the sub-pathways distribution as well
table(dat.met$`SUB PATHWAY`, useNA = "ifany")
#Check the number of sub-pathways
length(table(dat.met$`SUB PATHWAY`, useNA = "ifany")) #114
#check the range of metabolite numbers falling under different sub-pathways
range(table(dat.met$`SUB PATHWAY`, useNA = "ifany")) #(1, 336)


#Loading nearest neighbors (K=10) imputed COPDGene Data -- 20Aug2018
load("D:/OneDrive - The University of Colorado Denver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/OneDrive_2021-05-19/MetaboGuru_Data/COPDGene_Metabolon/copd_metabolon_imputed_20Aug2018.Rda")
#store the data as a data-frame - subjects in rows, metabolites in columns
dat.sub = as.data.frame(copd_metabolon)
dim(dat.sub) #1298 x 1035
colnames(dat.sub)

#select subjects/patients only from the Visit-2
dat.sub.v2 = as.data.frame(dat.sub[which(dat.sub$TIME_POINT == "Visit 2"), ])
#check number of subjects/patients from Visit-2
dim(dat.sub.v2) #1136 x 1035

COPDGene.subset = read.table("D:/OneDrive - The University of Colorado Denver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/OneDrive_2021-05-19/MetaboGuru_Data/COPDGene_Metabolon/COPDGene_subset.txt",
                             sep = "", header = T)
dim(COPDGene.subset) #17129 x 14
colnames(COPDGene.subset)

COPDGene.subset.v2 = COPDGene.subset[which(COPDGene.subset$visitnum == 2),]
dim(COPDGene.subset.v2) #6758 x 14

#match with SIDs from the two data-frames
sid.overlap = which(COPDGene.subset.v2$sid %in% dat.sub.v2$SID)
length(sid.overlap)

#select the sample metadata from Visit-2
COPDGene.subset.v2.metadata = COPDGene.subset.v2[sid.overlap,]
dim(COPDGene.subset.v2.metadata) #1136 x 14
#pick the columns for SID, Visit #, age-at-visit, gender
clnames = which(colnames(COPDGene.subset.v2.metadata)  %in% c("sid", "visitnum", "age_visit", "gender"))
#check column names of the 4 picked columns
colnames(COPDGene.subset.v2.metadata)[clnames]
#change column-names to - SID, TIME_POINT, AGEE, GENDER - to match with those in metabolites metadata
colnames(COPDGene.subset.v2.metadata)[clnames] = c("SID", "TIME_POINT", "AGE", "GENDER")

#match SIDs in dat.sub.2 with those in sample-metadata
ind = match(dat.sub.v2$SID, COPDGene.subset.v2.metadata$SID)
#Check for if the above SID match is accurately done
samp = sample(nrow(dat.sub.v2), 10, replace = F)
#Check if "Age" matches
sum(dat.sub.v2$AGE[samp] != COPDGene.subset.v2.metadata$AGE[ind[samp]])
#Check if "Gender" matches
sum(dat.sub.v2$GENDER[samp] != COPDGene.subset.v2.metadata$GENDER[ind[samp]])
#Check if "BMI" matches
sum(dat.sub.v2$BMI[samp] != COPDGene.subset.v2.metadata$BMI[ind[samp]])
#Check if "ATS_PackYears" matches
sum(dat.sub.v2$ATS_PackYears[samp] != COPDGene.subset.v2.metadata$ATS_PackYears[ind[samp]],
    na.rm = T)

#Check dimensions and column-names of the selected metabolites' and samples' metadata
dim(dat.sub.v2) #1136 x 1035
dim(COPDGene.subset.v2.metadata) #1136 x 14
colnames(dat.sub.v2)
colnames(COPDGene.subset.v2.metadata)

#Create the final data-frame: merge the two selected data-frames
COPDGene.subset.v2.final = merge(dat.sub.v2, COPDGene.subset.v2.metadata, by = c("SID","BMI", "AGE", "GENDER", "ATS_PackYears"))
#check the dimension of the final data-frame
dim(COPDGene.subset.v2.final) #1136 x 1044
#check the first few column-names
head(colnames(COPDGene.subset.v2.final))

#Two outcomes considered for analysis - categorical (Gender) and continuous (FEV-1/FVC ratio)
#Outcome 1: Gender
out_gender = COPDGene.subset.v2.final$GENDER
table(out_gender) #1: 572, 2: 564
sum(is.na(out_gender)) #0

#Outcome 2: FEV1/FVC ratio
fevrat = COPDGene.subset.v2.final$FEV1_FVC_utah
#plot the histogram
hist(fevrat)
#check for NA's
nas_fevrat = which(is.na(fevrat))
#check number of NA's
length(nas_fevrat) #16
#convert FEV-1/FVC ratios to their corresponding logit transform
out_logFEVrat = log(fevrat/(1-fevrat))
#check the histogram of logit of FEV-1/FVC ratios
hist(out_logFEVrat, breaks = 15)

#Collect the "confounders" (covariates) for Outcome 2
conf_ind = which(colnames(COPDGene.subset.v2.final) %in%
                   c("AGE","GENDER","BMI","smoking_status","ATS_PackYears", "ccenter"))
conf_logFEVrat = COPDGene.subset.v2.final[, conf_ind]
#check dimension of the covariates data-frame
dim(conf_logFEVrat) #1136 x 6
#check the column names of the covariates data-frame
colnames(conf_logFEVrat)

#Check for NA's in the subjects for each explanatory variable
sumnas_pred = apply(conf_logFEVrat, 2, function(s) sum(is.na(s)))
#Check for NA's in the explanatory variables for each subject
sumnas_samp = apply(conf_logFEVrat, 1, function(s) sum(is.na(s)))

#check how many explanatory variables have at least one NA
sum(sumnas_pred > 0) #1
#which explanatory variable has NA's
colnames(conf_logFEVrat)[sumnas_pred > 0] #"ATS_PackYears"

#check how many subjects have at least one NA
sum(sumnas_samp > 0) #66
#which subjects (rows) have NA's
nas_samp = which(sumnas_samp > 0)
#display the coavariates info for these subjects
conf_logFEVrat[nas_samp, ]
#how many are such subjects
length(nas_samp) #66


#Omit the samples with NA's in either ATS_PackYears OR FEV1/FVC-ratio
#collect the row-id's of such subjects first
nas_fevrat_samp = union(nas_fevrat, nas_samp)
#how many are there
length(nas_fevrat_samp) #82
#Omit these subjects from the final data-frame
COPDGene.subset.v2.final2 = COPDGene.subset.v2.final[-nas_fevrat_samp,]
#check dimension of data-frame after the above subjects' omission
dim(COPDGene.subset.v2.final2) #1054 x 1044


#Final data-frames
#Covariates' data-frame - remove the subjects with NA's
conf_logFEVrat.final2 = conf_logFEVrat[-nas_fevrat_samp,]
#check dimension
dim(conf_logFEVrat.final2) #1054 x 6
#check column names
colnames(conf_logFEVrat.final2)

#Dichotomize the covariate - data collection center
#NJC: 0, UIA: 1
conf_logFEVrat.final2$ccenter = ifelse(conf_logFEVrat.final2$ccenter == "NJC", 0, 1)
#Dichotomize the covariate - gender
#Female: 0, Male: 1
conf_logFEVrat.final2$GENDER = ifelse(conf_logFEVrat.final2$GENDER == 2, 0, 1)
#Dichotomize the covariate - smoking-status
#Former: 0, Current: 1
conf_logFEVrat.final2$smoking_status = ifelse(conf_logFEVrat.final2$smoking_status == 1, 0, 1) 

#Outcomes - remove the subjects with NA's
#outcome - FEV-1/FVC ratio
fevrat.final2 = fevrat[-nas_fevrat_samp]
#logit of ratio
out_logFEVrat.final2 = out_logFEVrat[-nas_fevrat_samp]

#outcome - gender
out_gender.final2 = out_gender[-nas_fevrat_samp]

#Finalize the metabolites 
#check columns-names
colnames(COPDGene.subset.v2.final2)
#collect only the metabolites abundances
mets.all.final2 = COPDGene.subset.v2.final2[, 31:1035]
#take log-transformation of metabolite abundances
log.mets.all.final2 = log(mets.all.final2)
#check dimension of data-frame
dim(mets.all.final2) #1054 x 1005

#W=write metabolites data into a CSV-file
filename = "mets_p1005n1054_final2.csv"
write.csv(mets.all.final2, filename, 
          row.names = F)

#group metabolites as per their super-pathways' memberships
#match metabolite names with those within metabolite metadata
met_sub_super_path = match(colnames(COPDGene.subset.v2.final2)[31:1035], dat.met$BIOCHEMICAL) 
#check the number of metabolites
length(met_sub_super_path) #1005
#check the unique names of the super-pathways
sup.path.names = unique(dat.met$`SUPER PATHWAY`[met_sub_super_path])
#tabulate the metabolites in super-pathways groups
table(dat.met$`SUPER PATHWAY`[met_sub_super_path], useNA = "ifany")
#check number of super-pathways
nos.superPath = length(sup.path.names) #10

#collect indices of metabolites belonging to respective super-pathways in a list
sup.path.ind = vector(mode = "list", length = nos.superPath)
#names the elements of the list as per the super-pathway names 
names(sup.path.ind) = sup.path.names
#store the super-pathway grouping indices in respective list elements
for(k in 1:nos.superPath){
  sup.path.ind[[k]] = which(dat.met$`SUPER PATHWAY`[met_sub_super_path] %in% sup.path.names[k])
}

#Gather more info on the "Partially Characterized Molecules"
dat.met.1 = dat.met[met_sub_super_path,]
dim(dat.met.1) #1005 x 13
names(sup.path.ind)[9] #"Partially Characterized Molecules"

#names of the metabolites that are  partially characterized
dat.met.1$BIOCHEMICAL[unlist(sup.path.ind[[9]])] #"glycine conjugate of C10H12O2*"     "glycine conjugate of C10H14O2 (1)*"
#Compund IDs of the metabolites that are  partially characterized
dat.met.1$`COMP ID`[unlist(sup.path.ind[[9]])] #62145, 62146
#Chemical IDs of the metabolites that are  partially characterized
dat.met.1$`CHEMICAL ID`[unlist(sup.path.ind[[9]])] #100020253, 100020254

#cross-check once more the total number of metabolites
sum(unlist(lapply(sup.path.ind, length))) #1005

#Delete "Partially Characterized Molecules" and "NA" super-pathway categories
#define a new list for storing the final super-pathways
sup.path.ind.v2 = list()
#store the first 8 super-pathways
sup.path.ind.v2 = sup.path.ind[1:8]
#rename the super-pathways with shorters legends
names(sup.path.ind.v2) = c("Lipid","Xeno","Amino","Carb","CofacVit",
                           "Nucle","Energy", "Pep")

#create a new list to store indices of metabolites as per the finally selected super-pathways
mets.classes.final2 = vector(mode = "list", length = length(sup.path.ind.v2))
ind.v2.sub = c()
#Create data-frames of metabolites belonging to the selected 8 super-pathways
for(k in 1:length(sup.path.ind.v2)){
  class.ind = which(colnames(mets.all.final2)  %in% dat.met$BIOCHEMICAL[met_sub_super_path][unlist(sup.path.ind.v2[k])])
  print(sum(class.ind != unlist(sup.path.ind.v2[k])))
  mets.classes.final2[[k]] = mets.all.final2[, class.ind]
  filename = paste(c(names(sup.path.ind.v2)[k], ".csv"), collapse = "")
  write.csv(mets.all.final2[, class.ind], filename, 
            row.names = F)
  ind.v2.sub = c(ind.v2.sub, unlist(sup.path.ind.v2[[k]]))
}

#8 super-pathways grouped side-by-side
names(sup.path.ind.v2) #"Lipid"    "Xeno"     "Amino"    "Carb"     "CofacVit" "Nucle"    "Energy"   "Pep"  
#final selected metabolites data-frame with 8 super-pathways grouped side-by-side 
mets.p761.final2 = mets.all.final2[, ind.v2.sub]
#check dimension
dim(mets.p761.final2) #1054 x 761

#write the final metabolites data-frame in a CSV-file, used in the analysis later
filename = "mets_allsupPath_final2_p761.csv"
write.csv(mets.p761.final2, filename, 
          row.names = F)


save.image("dataProcessing_COPDGene_Meff_Final_061324.RData")





