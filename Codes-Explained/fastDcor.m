%fastDcor, based on the function fastDcov
%A fast algorithm for computing distance correlation
%Authors: Arin Chauduri & Wenhao Hu
%Link to teh research Gate article: 
%https://www.researchgate.net/publication/328575900_A_fast_algorithm_for_computing_distance_correlation#fullTextFileContent

function dcor = fastDcor(x,y)
dcovxy = fastDcov(x,y);
dcovx = fastDcov(x,x);
dcovy = fastDcov(y,y);

dcor2 = dcovxy/sqrt(dcovx * dcovy);

dcor = sqrt(dcor2);
end



