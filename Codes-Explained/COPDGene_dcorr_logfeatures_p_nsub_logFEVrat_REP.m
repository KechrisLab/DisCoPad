%MWSL
%Retained p = M features (metabolites) from the total set of 1005 metabolites; "NA" and
%"Partially Characterized Molecules" omitted

%store the final metabolites dataset prepared in the data processing step
file1 = 'mets_allsupPath_final2_p761.csv';

%read the data in the 'features' variable
features = readmatrix(file1);
%check dimension
size(features)

%Compute distance correlation for REP subsets of p features
%get the sample / subject indices stored for Pearson correlation
%computation
file2 = 'sampSel_nsub-p_REP.csv';

 
sampSel = readmatrix(file2);
%metsRand = readmatrix(file3);

rep=size(sampSel,1);
p_sub = size(features,2);

%Use parallel for loop (parfor) to compute the distance correlation
%Only compute the upper traingular matrix

tic
parfor r = 1:rep
   r
   d_corr = ones(p_sub,p_sub);
   file3 = strcat('logfeat_sub_res_nsub-p_', num2str(r),'.csv');
   logfeatures_sub = readmatrix(file3);
   for k = 1:(p_sub-1)
    for l = (k+1):p_sub
        %[k, l]
        d_corr(k,l) = fastDcor(logfeatures_sub(:,k), logfeatures_sub(:,l));
        d_corr(l,k) = d_corr(k,l);
    end
   end 
   filename = strcat('dcorr_logfeatures_sub_RES_p-nsub_logFEVrat_', num2str(r));
   writematrix(d_corr, filename);
end
toc


%Check for positive definite
symCheck=zeros(rep, 1);
nearPD_dcorr_logfeatures_sub=ones(p_sub,p_sub,rep);

for k = 1:rep
    k
    filename = strcat('dcorr_logfeatures_sub_RES_p-nsub_logFEVrat_', num2str(k));
    dcorr_logfeatures_sub = readmatrix(filename);
    symCheck(k) = issymmetric(dcorr_logfeatures_sub);
    eVal = eig(dcorr_logfeatures_sub);

    %convert the DisCo matrix to its nearest SPD matrix
    if sum(eVal < 0) > 0
     nearPD_dcorr_logfeatures_sub(:,:,k) = nearestSPD(dcorr_logfeatures_sub); 
    else
     nearPD_dcorr_logfeatures_sub(:,:,k) = dcorr_logfeatures_sub;   
    end        
end
sum(symCheck) %should equal rep

pdCheck = zeros(rep,1);
symCheck2 = zeros(rep,1);

for l=1:rep
    l
    ch = chol(nearPD_dcorr_logfeatures_sub(:,:,l));
    pdCheck(l) = all(diag(ch) > 0);
    symCheck2(l) = issymmetric(nearPD_dcorr_logfeatures_sub(:,:,l));
end
sum(pdCheck) %should equal rep
sum(symCheck2) %should equal rep

for k=1:rep
    k
    filename = strcat('COPDGene_nearPD_dcorr_logfeatures_sub_RES_p-nsub_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_sub(:,:,k), filename);
end



%%%%% Compute distance correlations separately within each of the 8 metabolic pathways
%%%%%%%%%% Outcome: log-FEV-1/FVC-ratio
%p_sub = [365, 99, 181, 24, 25, 32, 10, 25];
rep = 100;

tic
parfor r = 1:rep
   
   r
    
   p_sub = [365, 99, 181, 24, 25, 32, 10, 25]; %cardinality of different metabolic super-pathways
    
   file1 = strcat("/Users/dnandy/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/R-Codes/n100p761_rep100_logFEVrat/logfeat_sub_res_n100p365class1_", num2str(r), ".csv");
   file2 = strcat("/Users/dnandy/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/R-Codes/n100p761_rep100_logFEVrat/logfeat_sub_res_n100p99class2_", num2str(r), ".csv");
   file3 = strcat("/Users/dnandy/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/R-Codes/n100p761_rep100_logFEVrat/logfeat_sub_res_n100p181class3_", num2str(r), ".csv");
   file4 = strcat("/Users/dnandy/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/R-Codes/n100p761_rep100_logFEVrat/logfeat_sub_res_n100p24class4_", num2str(r), ".csv");
   file5 = strcat("/Users/dnandy/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/R-Codes/n100p761_rep100_logFEVrat/logfeat_sub_res_n100p25class5_", num2str(r), ".csv");
   file6 = strcat("/Users/dnandy/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/R-Codes/n100p761_rep100_logFEVrat/logfeat_sub_res_n100p32class6_", num2str(r), ".csv");
   file7 = strcat("/Users/dnandy/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/R-Codes/n100p761_rep100_logFEVrat/logfeat_sub_res_n100p10class7_", num2str(r), ".csv");
   file8 = strcat("/Users/dnandy/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Postdoc_Project/Effective-#-Of-Tests-Metabolomics/R-Codes/n100p761_rep100_logFEVrat/logfeat_sub_res_n100p25class8_", num2str(r), ".csv");

   mat1 = readmatrix(file1);
   mat2 = readmatrix(file2);
   mat3 = readmatrix(file3);
   mat4 = readmatrix(file4);
   mat5 = readmatrix(file5);
   mat6 = readmatrix(file6);
   mat7 = readmatrix(file7);
   mat8 = readmatrix(file8);
   
   %Class 1
   [r,1]
   dcorr_logfeatures_class_sub = ones(p_sub(1));
   for k = 1:(p_sub(1)-1)
       %k
    for l = (k+1):p_sub(1)
        %[k, l]
        dcorr_logfeatures_class_sub(k,l) = fastDcor(mat1(:,k), mat1(:,l));
        dcorr_logfeatures_class_sub(l,k) = dcorr_logfeatures_class_sub(k,l);
    end
   end 
   filename = strcat('dcorr_logfeatures_class1_sub_', num2str(r));
   writematrix(dcorr_logfeatures_class_sub, filename);

   %Class 2
   [r,2]
   dcorr_logfeatures_class_sub = ones(p_sub(2));
   for k = 1:(p_sub(2)-1)
       %k
    for l = (k+1):p_sub(2)
        %[k, l]
        dcorr_logfeatures_class_sub(k,l) = fastDcor(mat2(:,k), mat2(:,l));
        dcorr_logfeatures_class_sub(l,k) = dcorr_logfeatures_class_sub(k,l);
    end
   end 
   filename = strcat('dcorr_logfeatures_class2_sub_', num2str(r));
   writematrix(dcorr_logfeatures_class_sub, filename);
   
   %Class 3
   [r,3]
   dcorr_logfeatures_class_sub = ones(p_sub(3));
   for k = 1:(p_sub(3)-1)
       %k
    for l = (k+1):p_sub(3)
        %[k, l]
        dcorr_logfeatures_class_sub(k,l) = fastDcor(mat3(:,k), mat3(:,l));
        dcorr_logfeatures_class_sub(l,k) = dcorr_logfeatures_class_sub(k,l);
    end
   end 
   filename = strcat('dcorr_logfeatures_class3_sub_', num2str(r));
   writematrix(dcorr_logfeatures_class_sub, filename);
   
   %Class 4
   [r,4]
   dcorr_logfeatures_class_sub = ones(p_sub(4));
   for k = 1:(p_sub(4)-1)
       %k
    for l = (k+1):p_sub(4)
        %[k, l]
        dcorr_logfeatures_class_sub(k,l) = fastDcor(mat4(:,k), mat4(:,l));
        dcorr_logfeatures_class_sub(l,k) = dcorr_logfeatures_class_sub(k,l);
    end
   end 
   filename = strcat('dcorr_logfeatures_class4_sub_', num2str(r));
   writematrix(dcorr_logfeatures_class_sub, filename);
   
   %Class 5
   [r,5]
   dcorr_logfeatures_class_sub = ones(p_sub(5));
   for k = 1:(p_sub(5)-1)
       %k
    for l = (k+1):p_sub(5)
        %[k, l]
        dcorr_logfeatures_class_sub(k,l) = fastDcor(mat5(:,k), mat5(:,l));
        dcorr_logfeatures_class_sub(l,k) = dcorr_logfeatures_class_sub(k,l);
    end
   end 
   filename = strcat('dcorr_logfeatures_class5_sub_', num2str(r));
   writematrix(dcorr_logfeatures_class_sub, filename);
   
   %Class 6
   [r,6]
   dcorr_logfeatures_class_sub = ones(p_sub(6));
   for k = 1:(p_sub(6)-1)
       %k
    for l = (k+1):p_sub(6)
        %[k, l]
        dcorr_logfeatures_class_sub(k,l) = fastDcor(mat6(:,k), mat6(:,l));
        dcorr_logfeatures_class_sub(l,k) = dcorr_logfeatures_class_sub(k,l);
    end
   end
   filename = strcat('dcorr_logfeatures_class6_sub_', num2str(r));
   writematrix(dcorr_logfeatures_class_sub, filename);
   
   %Class 7
   [r,7]
   dcorr_logfeatures_class_sub = ones(p_sub(7));
   for k = 1:(p_sub(7)-1)
       %k
    for l = (k+1):p_sub(7)
        %[k, l]
        dcorr_logfeatures_class_sub(k,l) = fastDcor(mat7(:,k), mat7(:,l));
        dcorr_logfeatures_class_sub(l,k) = dcorr_logfeatures_class_sub(k,l);
    end
   end 
   filename = strcat('dcorr_logfeatures_class7_sub_', num2str(r));
   writematrix(dcorr_logfeatures_class_sub, filename);
   
    %Class 8
    [r,8]
    dcorr_logfeatures_class_sub = ones(p_sub(8));
   for k = 1:(p_sub(8)-1)
       %k
    for l = (k+1):p_sub(8)
        %[k, l]
        dcorr_logfeatures_class_sub(k,l) = fastDcor(mat8(:,k), mat8(:,l));
        dcorr_logfeatures_class_sub(l,k) = dcorr_logfeatures_class_sub(k,l);
    end
   end 
   filename = strcat('dcorr_logfeatures_class8_sub_', num2str(r));
   writematrix(dcorr_logfeatures_class_sub, filename);

end
toc

p_sub = [365, 99, 181, 24, 25, 32, 10, 25];
%Check
symCheck=zeros(rep, 8);
nearPD_dcorr_logfeatures_class1_sub=ones(p_sub(1),p_sub(1),rep);
nearPD_dcorr_logfeatures_class2_sub=ones(p_sub(2),p_sub(2),rep);
nearPD_dcorr_logfeatures_class3_sub=ones(p_sub(3),p_sub(3),rep);
nearPD_dcorr_logfeatures_class4_sub=ones(p_sub(4),p_sub(4),rep);
nearPD_dcorr_logfeatures_class5_sub=ones(p_sub(5),p_sub(5),rep);
nearPD_dcorr_logfeatures_class6_sub=ones(p_sub(6),p_sub(6),rep);
nearPD_dcorr_logfeatures_class7_sub=ones(p_sub(7),p_sub(7),rep);
nearPD_dcorr_logfeatures_class8_sub=ones(p_sub(8),p_sub(8),rep);

for k = 1:rep
    k
    filename = strcat('dcorr_logfeatures_class1_sub_', num2str(k));
    dcorr_logfeatures_class1_sub = readmatrix(filename);
    symCheck(k,1) = issymmetric(dcorr_logfeatures_class1_sub);

    filename = strcat('dcorr_logfeatures_class2_sub_', num2str(k));
    dcorr_logfeatures_class2_sub = readmatrix(filename);
    symCheck(k,2) = issymmetric(dcorr_logfeatures_class2_sub);

    filename = strcat('dcorr_logfeatures_class3_sub_', num2str(k));
    dcorr_logfeatures_class3_sub = readmatrix(filename);
    symCheck(k,3) = issymmetric(dcorr_logfeatures_class3_sub);

    filename = strcat('dcorr_logfeatures_class4_sub_', num2str(k));
    dcorr_logfeatures_class4_sub = readmatrix(filename);
    symCheck(k,4) = issymmetric(dcorr_logfeatures_class4_sub);

    filename = strcat('dcorr_logfeatures_class5_sub_', num2str(k));
    dcorr_logfeatures_class5_sub = readmatrix(filename);
    symCheck(k,5) = issymmetric(dcorr_logfeatures_class5_sub);

    filename = strcat('dcorr_logfeatures_class6_sub_', num2str(k));
    dcorr_logfeatures_class6_sub = readmatrix(filename);
    symCheck(k,6) = issymmetric(dcorr_logfeatures_class6_sub);

    filename = strcat('dcorr_logfeatures_class7_sub_', num2str(k));
    dcorr_logfeatures_class7_sub = readmatrix(filename);
    symCheck(k,7) = issymmetric(dcorr_logfeatures_class7_sub);

    filename = strcat('dcorr_logfeatures_class8_sub_', num2str(k));
    dcorr_logfeatures_class8_sub = readmatrix(filename);
    symCheck(k,8) = issymmetric(dcorr_logfeatures_class8_sub);
    
    eVal1 = eig(dcorr_logfeatures_class1_sub);
    eVal2 = eig(dcorr_logfeatures_class2_sub);
    eVal3 = eig(dcorr_logfeatures_class3_sub);
    eVal4 = eig(dcorr_logfeatures_class4_sub);
    eVal5 = eig(dcorr_logfeatures_class5_sub);
    eVal6 = eig(dcorr_logfeatures_class6_sub);
    eVal7 = eig(dcorr_logfeatures_class7_sub);
    eVal8 = eig(dcorr_logfeatures_class8_sub);
    
    %convert the super-pathways grouped DisCo matrices to their nearest SPD
    %matrices

        if sum(eVal1 < 0) > 0
         nearPD_dcorr_logfeatures_class1_sub(:,:,k) = nearestSPD(dcorr_logfeatures_class1_sub); 
        else
         nearPD_dcorr_logfeatures_class1_sub(:,:,k) = dcorr_logfeatures_class1_sub;  
        end        
        
        if sum(eVal2 < 0) > 0
     nearPD_dcorr_logfeatures_class2_sub(:,:,k) = nearestSPD(dcorr_logfeatures_class2_sub);
        else
     nearPD_dcorr_logfeatures_class2_sub(:,:,k) = dcorr_logfeatures_class2_sub;  
        end 
       
        if sum(eVal3 < 0) > 0
     nearPD_dcorr_logfeatures_class3_sub(:,:,k) = nearestSPD(dcorr_logfeatures_class3_sub); 
        else
     nearPD_dcorr_logfeatures_class3_sub(:,:,k) = dcorr_logfeatures_class3_sub;  
        end 

        if sum(eVal4 < 0) > 0
     nearPD_dcorr_logfeatures_class4_sub(:,:,k) = nearestSPD(dcorr_logfeatures_class4_sub); 
        else
     nearPD_dcorr_logfeatures_class4_sub(:,:,k) = dcorr_logfeatures_class4_sub;  
        end 

        if sum(eVal5 < 0) > 0
     nearPD_dcorr_logfeatures_class5_sub(:,:,k) = nearestSPD(dcorr_logfeatures_class5_sub); 
        else
     nearPD_dcorr_logfeatures_class5_sub(:,:,k) = dcorr_logfeatures_class5_sub;  
        end 

        if sum(eVal6 < 0) > 0
     nearPD_dcorr_logfeatures_class6_sub(:,:,k) = nearestSPD(dcorr_logfeatures_class6_sub); 
        else
     nearPD_dcorr_logfeatures_class6_sub(:,:,k) = dcorr_logfeatures_class6_sub;  
        end

        if sum(eVal7 < 0) > 0
     nearPD_dcorr_logfeatures_class7_sub(:,:,k) = nearestSPD(dcorr_logfeatures_class7_sub); 
        else
     nearPD_dcorr_logfeatures_class7_sub(:,:,k) = dcorr_logfeatures_class7_sub;  
        end

        if sum(eVal8 < 0) > 0
     nearPD_dcorr_logfeatures_class8_sub(:,:,k) = nearestSPD(dcorr_logfeatures_class8_sub); 
        else
     nearPD_dcorr_logfeatures_class8_sub(:,:,k) = dcorr_logfeatures_class8_sub;  
        end 
      
end
sum(symCheck) %should equal REP

%check for PD
pdCheck = zeros(rep,8);
symCheck2 = zeros(rep,8);
for l=1:rep
    l
    ch1 = chol(nearPD_dcorr_logfeatures_class1_sub(:,:,l));
    pdCheck(l,1) = all(diag(ch1) > 0);
    symCheck2(l,1) = issymmetric(nearPD_dcorr_logfeatures_class1_sub(:,:,l));
    
    ch2 = chol(nearPD_dcorr_logfeatures_class2_sub(:,:,l));
    pdCheck(l,2) = all(diag(ch2) > 0);
    symCheck2(l,2) = issymmetric(nearPD_dcorr_logfeatures_class2_sub(:,:,l));
    
    ch3 = chol(nearPD_dcorr_logfeatures_class3_sub(:,:,l));
    pdCheck(l,3) = all(diag(ch3) > 0);
    symCheck2(l,3) = issymmetric(nearPD_dcorr_logfeatures_class3_sub(:,:,l));
    
    ch4 = chol(nearPD_dcorr_logfeatures_class4_sub(:,:,l));
    pdCheck(l,4) = all(diag(ch4) > 0);
    symCheck2(l,4) = issymmetric(nearPD_dcorr_logfeatures_class4_sub(:,:,l));
    
    ch5 = chol(nearPD_dcorr_logfeatures_class5_sub(:,:,l));
    pdCheck(l,5) = all(diag(ch5) > 0);
    symCheck2(l,5) = issymmetric(nearPD_dcorr_logfeatures_class5_sub(:,:,l));
    
    ch6 = chol(nearPD_dcorr_logfeatures_class6_sub(:,:,l));
    pdCheck(l,6) = all(diag(ch6) > 0);
    symCheck2(l,6) = issymmetric(nearPD_dcorr_logfeatures_class6_sub(:,:,l));
    
    ch7 = chol(nearPD_dcorr_logfeatures_class7_sub(:,:,l));
    pdCheck(l,7) = all(diag(ch7) > 0);
    symCheck2(l,7) = issymmetric(nearPD_dcorr_logfeatures_class7_sub(:,:,l));
    
    ch8 = chol(nearPD_dcorr_logfeatures_class8_sub(:,:,l));
    pdCheck(l,8) = all(diag(ch8) > 0);
    symCheck2(l,8) = issymmetric(nearPD_dcorr_logfeatures_class8_sub(:,:,l));
    
end
sum(pdCheck) %should equal rep
sum(symCheck2) %should equal rep


for k=1:rep

    k

    filename1 = strcat('COPDGene_nearPD_dcorr_logfeatures_p365nsub_class1_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_class1_sub(:,:,k), filename1);
    
    filename2 = strcat('COPDGene_nearPD_dcorr_logfeatures_p99nsub_class2_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_class2_sub(:,:,k), filename2);
    
    filename3 = strcat('COPDGene_nearPD_dcorr_logfeatures_p181nsub_class3_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_class3_sub(:,:,k), filename3);
    
    filename4 = strcat('COPDGene_nearPD_dcorr_logfeatures_p24nsub_class4_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_class4_sub(:,:,k), filename4);
    
    filename5 = strcat('COPDGene_nearPD_dcorr_logfeatures_p25nsub_class5_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_class5_sub(:,:,k), filename5);
    
    filename6 = strcat('COPDGene_nearPD_dcorr_logfeatures_p32nsub_class6_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_class6_sub(:,:,k), filename6);
    
    filename7 = strcat('COPDGene_nearPD_dcorr_logfeatures_p10nsub_class7_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_class7_sub(:,:,k), filename7);
    
    filename8 = strcat('COPDGene_nearPD_dcorr_logfeatures_p25nsub_class8_logFEVrat_', num2str(k));
    writematrix(nearPD_dcorr_logfeatures_class8_sub(:,:,k), filename8);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random grouping of metabolites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Random Groupings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randRep = 1; % ONE random grouping; Change randRep to any number to get that many ramdom grouping results

filename = 'mets_allsupPath_final2_p761.csv';
features = readmatrix(filename);
size(features)
p_sub = size(features,2);
p_sub_8class = [365, 99, 181, 24, 25, 32, 10, 25];
p_sub_8class_cumSum = cumsum(p_sub_8class);

rep = REP;
symCheck =  zeros(rep,8,randRep);
pdCheck  =  zeros(rep,8,randRep);
symCheck2 = zeros(rep,8, randRep);

%Random grouping just once
for rand=1:randRep
    
    s = RandStream('mlfg6331_64'); 
    p_sub_randInd = datasample(s, 1:p_sub, p_sub, 'Replace', false);
    writematrix(p_sub_randInd, 'rand0Ind_nsub-p_Final');  

    nearPD_dcorr_logfeatures_randclass1_sub=ones(p_sub_8class(1),p_sub_8class(1),rep);
    nearPD_dcorr_logfeatures_randclass2_sub=ones(p_sub_8class(2),p_sub_8class(2),rep);
    nearPD_dcorr_logfeatures_randclass3_sub=ones(p_sub_8class(3),p_sub_8class(3),rep);
    nearPD_dcorr_logfeatures_randclass4_sub=ones(p_sub_8class(4),p_sub_8class(4),rep);
    nearPD_dcorr_logfeatures_randclass5_sub=ones(p_sub_8class(5),p_sub_8class(5),rep);
    nearPD_dcorr_logfeatures_randclass6_sub=ones(p_sub_8class(6),p_sub_8class(6),rep);
    nearPD_dcorr_logfeatures_randclass7_sub=ones(p_sub_8class(7),p_sub_8class(7),rep);
    nearPD_dcorr_logfeatures_randclass8_sub=ones(p_sub_8class(8),p_sub_8class(8),rep);


    tic
    for r = 1:rep
   
       strcat('Rand= ', num2str(rand), ' | Rep= ', num2str(r))

       filename = strcat("dcorr_logfeatures_sub_res_p-nsub_logFEVrat_", num2str(r), ".txt");
       dcorrMatp761 = readmatrix(filename);

       %Random group 1
       dcorr_randclass1 = dcorrMatp761(p_sub_randInd(1:p_sub_8class_cumSum(1)), p_sub_randInd(1:p_sub_8class_cumSum(1)));
       symCheck(r,1,rand) = issymmetric(dcorr_randclass1);
       
       eVal1 = eig(dcorr_randclass1);    
       if sum(eVal1 < 0) > 0
           nearPD_dcorr_logfeatures_randclass1_sub(:,:,r) = nearestSPD(dcorr_randclass1); 
       else
           nearPD_dcorr_logfeatures_randclass1_sub(:,:,r) = dcorr_randclass1;  
       end
       
       ch = chol(nearPD_dcorr_logfeatures_randclass1_sub(:,:,r));
       pdCheck(r,1,rand) = all(diag(ch) > 0);
       symCheck2(r,1,rand) = issymmetric(nearPD_dcorr_logfeatures_randclass1_sub(:,:,r));

       %Random group 2
       dcorr_randclass2 = dcorrMatp761(p_sub_randInd(p_sub_8class_cumSum(1)+1:p_sub_8class_cumSum(2)), p_sub_randInd(p_sub_8class_cumSum(1)+1:p_sub_8class_cumSum(2)));
       symCheck(r,2,rand) = issymmetric(dcorr_randclass2);
       
       eVal2 = eig(dcorr_randclass2);    
       if sum(eVal2 < 0) > 0
           nearPD_dcorr_logfeatures_randclass2_sub(:,:,r) = nearestSPD(dcorr_randclass2); 
       else
           nearPD_dcorr_logfeatures_randclass2_sub(:,:,r) = dcorr_randclass2;  
       end

       ch = chol(nearPD_dcorr_logfeatures_randclass2_sub(:,:,r));
       pdCheck(r,2,rand) = all(diag(ch) > 0);
       symCheck2(r,2,rand) = issymmetric(nearPD_dcorr_logfeatures_randclass2_sub(:,:,r));

       %Random group 3
       dcorr_randclass3 = dcorrMatp761(p_sub_randInd(p_sub_8class_cumSum(2)+1:p_sub_8class_cumSum(3)), p_sub_randInd(p_sub_8class_cumSum(2)+1:p_sub_8class_cumSum(3)));
       symCheck(r,3,rand) = issymmetric(dcorr_randclass3);
       
       eVal3 = eig(dcorr_randclass3);    
       if sum(eVal3 < 0) > 0
           nearPD_dcorr_logfeatures_randclass3_sub(:,:,r) = nearestSPD(dcorr_randclass3); 
       else
           nearPD_dcorr_logfeatures_randclass3_sub(:,:,r) = dcorr_randclass3;  
       end

       ch = chol(nearPD_dcorr_logfeatures_randclass3_sub(:,:,r));
       pdCheck(r,3,rand) = all(diag(ch) > 0);
       symCheck2(r,3,rand) = issymmetric(nearPD_dcorr_logfeatures_randclass3_sub(:,:,r));


       %Random group 4
       dcorr_randclass4 = dcorrMatp761(p_sub_randInd(p_sub_8class_cumSum(3)+1:p_sub_8class_cumSum(4)), p_sub_randInd(p_sub_8class_cumSum(3)+1:p_sub_8class_cumSum(4)));
       symCheck(r,4,rand) = issymmetric(dcorr_randclass4);
       
       eVal4 = eig(dcorr_randclass4);    
       if sum(eVal1 < 0) > 0
           nearPD_dcorr_logfeatures_randclass4_sub(:,:,r) = nearestSPD(dcorr_randclass4); 
       else
           nearPD_dcorr_logfeatures_randclass4_sub(:,:,r) = dcorr_randclass4;  
       end

       ch = chol(nearPD_dcorr_logfeatures_randclass4_sub(:,:,r));
       pdCheck(r,4,rand) = all(diag(ch) > 0);
       symCheck2(r,4,rand) = issymmetric(nearPD_dcorr_logfeatures_randclass4_sub(:,:,r));

       %Random group 5
       dcorr_randclass5 = dcorrMatp761(p_sub_randInd(p_sub_8class_cumSum(4)+1:p_sub_8class_cumSum(5)), p_sub_randInd(p_sub_8class_cumSum(4)+1:p_sub_8class_cumSum(5)));
       symCheck(r,5,rand) = issymmetric(dcorr_randclass5);
       
       eVal5 = eig(dcorr_randclass5);    
       if sum(eVal5 < 0) > 0
           nearPD_dcorr_logfeatures_randclass5_sub(:,:,r) = nearestSPD(dcorr_randclass5); 
       else
           nearPD_dcorr_logfeatures_randclass5_sub(:,:,r) = dcorr_randclass5;  
       end

       ch = chol(nearPD_dcorr_logfeatures_randclass5_sub(:,:,r));
       pdCheck(r,5,rand) = all(diag(ch) > 0);
       symCheck2(r,5,rand) = issymmetric(nearPD_dcorr_logfeatures_randclass5_sub(:,:,r));

       %Random group 6
       dcorr_randclass6 = dcorrMatp761(p_sub_randInd(p_sub_8class_cumSum(5)+1:p_sub_8class_cumSum(6)), p_sub_randInd(p_sub_8class_cumSum(5)+1:p_sub_8class_cumSum(6)));
       symCheck(r,6,rand) = issymmetric(dcorr_randclass6);
       
       eVal6 = eig(dcorr_randclass6);    
       if sum(eVal6 < 0) > 0
           nearPD_dcorr_logfeatures_randclass6_sub(:,:,r) = nearestSPD(dcorr_randclass6); 
       else
           nearPD_dcorr_logfeatures_randclass6_sub(:,:,r) = dcorr_randclass6;  
       end

       ch = chol(nearPD_dcorr_logfeatures_randclass6_sub(:,:,r));
       pdCheck(r,6,rand) = all(diag(ch) > 0);
       symCheck2(r,6,rand) = issymmetric(nearPD_dcorr_logfeatures_randclass6_sub(:,:,r));


       %Random group 7
       dcorr_randclass7 = dcorrMatp761(p_sub_randInd(p_sub_8class_cumSum(6)+1:p_sub_8class_cumSum(7)), p_sub_randInd(p_sub_8class_cumSum(6)+1:p_sub_8class_cumSum(7)));
       symCheck(r,7,rand) = issymmetric(dcorr_randclass7);
       
       eVal7 = eig(dcorr_randclass7);    
       if sum(eVal7 < 0) > 0
           nearPD_dcorr_logfeatures_randclass7_sub(:,:,r) = nearestSPD(dcorr_randclass7); 
       else
           nearPD_dcorr_logfeatures_randclass7_sub(:,:,r) = dcorr_randclass7;  
       end

       ch = chol(nearPD_dcorr_logfeatures_randclass7_sub(:,:,r));
       pdCheck(r,7,rand) = all(diag(ch) > 0);
       symCheck2(r,7,rand) = issymmetric(nearPD_dcorr_logfeatures_randclass7_sub(:,:,r));


       %Random group 8
       dcorr_randclass8 = dcorrMatp761(p_sub_randInd(p_sub_8class_cumSum(7)+1:p_sub_8class_cumSum(8)), p_sub_randInd(p_sub_8class_cumSum(7)+1:p_sub_8class_cumSum(8)));
       symCheck(r,8,rand) = issymmetric(dcorr_randclass8);
       
       eVal8 = eig(dcorr_randclass8);    
       if sum(eVal8 < 0) > 0
           nearPD_dcorr_logfeatures_randclass8_sub(:,:,r) = nearestSPD(dcorr_randclass8); 
       else
           nearPD_dcorr_logfeatures_randclass8_sub(:,:,r) = dcorr_randclass8;  
       end

       ch = chol(nearPD_dcorr_logfeatures_randclass8_sub(:,:,r));
       pdCheck(r,8,rand) = all(diag(ch) > 0);
       symCheck2(r,8,rand) = issymmetric(nearPD_dcorr_logfeatures_randclass8_sub(:,:,r));


    filename = strcat('Random/COPDGene_nearPD_dcorr_logfeatures_p365nsub_randclass1_logFEVrat_rand', num2str(rand-1), '_', num2str(r));
    writematrix(nearPD_dcorr_logfeatures_randclass1_sub(:,:,r), filename);
    
    filename = strcat('Random/COPDGene_nearPD_dcorr_logfeatures_p99nsub_randclass2_logFEVrat_rand', num2str(rand-1), '_', num2str(r));
    writematrix(nearPD_dcorr_logfeatures_randclass2_sub(:,:,r), filename);
    
    filename = strcat('Random/COPDGene_nearPD_dcorr_logfeatures_p181nsub_randclass3_logFEVrat_rand', num2str(rand-1), '_',num2str(r));
    writematrix(nearPD_dcorr_logfeatures_randclass3_sub(:,:,r), filename);
    
    filename = strcat('Random/COPDGene_nearPD_dcorr_logfeatures_p24nsub_randclass4_logFEVrat_rand', num2str(rand-1), '_', num2str(r));
    writematrix(nearPD_dcorr_logfeatures_randclass4_sub(:,:,r), filename);
    
    filename = strcat('Random/COPDGene_nearPD_dcorr_logfeatures_p25nsub_randclass5_logFEVrat_rand', num2str(rand-1), '_', num2str(r));
    writematrix(nearPD_dcorr_logfeatures_randclass5_sub(:,:,r), filename);
    
    filename = strcat('Random/COPDGene_nearPD_dcorr_logfeatures_p32nsub_randclass6_logFEVrat_rand', num2str(rand-1), '_', num2str(r));
    writematrix(nearPD_dcorr_logfeatures_randclass6_sub(:,:,r), filename);
    
    filename = strcat('Random/COPDGene_nearPD_dcorr_logfeatures_p10nsub_randclass7_logFEVrat_rand', num2str(rand-1), '_', num2str(r));
    writematrix(nearPD_dcorr_logfeatures_randclass7_sub(:,:,r), filename);
    
    filename = strcat('Random/COPDGene_nearPD_dcorr_logfeatures_p25nsub_randclass8_logFEVrat_rand', num2str(rand-1), '_', num2str(r));
    writematrix(nearPD_dcorr_logfeatures_randclass8_sub(:,:,r), filename);
   end
end


save COPDGene_dcorr_logfeatures_p_nsub_logFEVrat_REP.mat

