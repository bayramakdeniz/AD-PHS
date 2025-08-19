
% to activate glmnet
cd('/home/eadbwp1/shahram_bayram/ofrei/glmnet-matlab')
mex -setup fortran
addpath('/home/eadbwp1/shahram_bayram/ofrei/glmnet-matlab')
addpath('/home/eadbwp1/shahram_bayram/For_test')
addpath('/home/eadbwp1/shahram_bayram/github_phs_ofrei/EADB_PHS_2023')
cd('/home/eadbwp1/shahram_bayram/eadb_gsa_core_topmed/dbsnp156')
cd('/home/eadbwp1/shahram_bayram/For_test/AD-PHS-main/demo_example')


% generate new PHS model from genotype / age data:
load new_PHS_input_example % some junk data I made for testing

options = glmnetSet;
options.dfmax=500;
options.intr=false;  % silence the warning about intercept; it's not needed for coxph fit

% http://cran.nexr.com/web/packages/glmnet/vignettes/Coxnet.pdf
% y is an n length vector of failure/censoring times, and status is an n
% length vector with each entry, a 1 or a 0, indicating whether the corresponding
% entry in y is indicative of a failure time or right censoring time (1 for failure, 0
% for censoring)
y=cat(2, q.age, q.casevec);
nfolds=5;

nanind=find((isnan(y(:,1)) | isnan(y(:,2)))==0);
y1=y(nanind,:);
G=q.gmat(nanind,:);

    fit1=cvglmnet(G, y1, 'cox', options, 'deviance', nfolds);
    [~, indlam]=min( abs(fit1.lambda - fit1.lambda_min) );
    betas=fit1.glmnet_fit.beta(:, indlam); % effect size of the variants

