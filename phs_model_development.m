% prerequisites and environment on ZEUS server (EADB)
% wget https://hastie.su.domains/glmnet_matlab/glmnet_matlab_new.zip && unzip glmnet_matlab_new.zip
% ssh interactive
% module load matlab/R2019b

% options
alpha=0.03;  % based on cross-validation in EADB, alpha=0.03 seem to be a good generic choice for ElasticNet setting (LASSO vs RIDGE trade-off) 
nfolds = 5;
age_thresh = 60;
dfmax = 1500;  % max number of SNPs to consider while fitting GLMNET

% file inputs
clinical_file = "/home/eadbwp1/shahram_bayram/clinical_data_ofrei_2024_nov_12.tsv";
%genome_input='/home/eadbwp1/shahram_bayram/eadb_gsa_core_topmed/dbsnp156_withGERAD/eadb_core_pval_1e-4_common3snps';
genome_input='/home/eadbwp1/shahram_bayram/For_test/eadb_core_pval_1e-5';
% outputs:
%out_prefix='/home/eadbwp1/shahram_bayram/github_phs_ofrei/EADB_PHS_2023/eadb_ALZ_phs_apoe_genome_2025_01_30_v1';
out_prefix='/home/eadbwp1/shahram_bayram/For_test/of_bagen';
% history of previous runs:
% eadb_ALZ_phs_apoe_genome_2025_01_15_v1 - consider all EADB SNPs at pval=1e-5 (M=7246 SNPs)
% eadb_ALZ_phs_apoe_genome_2025_01_15_v2 - consider EADB SNPs overlapping with MVP and ADNI at pval=1e-5 (M=6956 SNPs)
% eadb_ALZ_phs_apoe_genome_2025_01_30_v1 - EADB SNPs overlapping with MVP and ADNI at pval=1e-4, nfold=10, age_thresh=55, not filter on in_GERAD; fitted from 18648 cases and 17783 controls

cd('/home/eadbwp1/shahram_bayram/ofrei/glmnet-matlab')
mex -setup fortran
addpath('/home/eadbwp1/shahram_bayram/ofrei/glmnet-matlab')
addpath('/home/eadbwp1/shahram_bayram/For_test')
addpath('/home/eadbwp1/shahram_bayram/github_phs_ofrei/EADB_PHS_2023')
cd('/home/eadbwp1/shahram_bayram/eadb_gsa_core_topmed/dbsnp156')
cd('/home/eadbwp1/shahram_bayram/For_test/')

% read clinical file
opts = detectImportOptions(clinical_file, "FileType", "text");
opts = setvartype(opts,{'IID'},'string');
opts = setvartype(opts,{'apoe'},'double');
t = readtable(clinical_file, opts); % 42140 subjects

% this age filtering (using only age_at_onset for cases) may work better after all, just to remind myself
t.('age')(:) = nan;

case_idx = t.('ad_status') == 1; 
t.('age')(case_idx) = str2double(t.('age_at_onset')(case_idx));
t.('age')(~case_idx) = nanmax(t.('age_last_exam')(~case_idx), str2double(t.('age_baseline')(~case_idx)));


t = t(isfinite(t.('age')), :);
t = t((t.('age') >= age_thresh), :);
t = t(isfinite(t.('apoe')), :);

time_var = t.('age');
event_var = t.('ad_status');  % 1 = cases, 0 = controls
id_var = t.('IID');
country_var = t.('country');
apoe_status = t.('apoe');

assert(all(isfinite(apoe_status)));
apoe_status2=zeros(length(apoe_status),3);
apoe_status2(find((apoe_status)==22),:)=[2 0 0].*ones(3,length(find((apoe_status)==22)))';
apoe_status2(find((apoe_status)==23),:)=[1 1 0].*ones(3,length(find((apoe_status)==23)))';
apoe_status2(find((apoe_status)==24),:)=[1 0 1].*ones(3,length(find((apoe_status)==24)))';
apoe_status2(find((apoe_status)==33),:)=[0 2 0].*ones(3,length(find((apoe_status)==33)))';
apoe_status2(find((apoe_status)==34),:)=[0 1 1].*ones(3,length(find((apoe_status)==34)))';
apoe_status2(find((apoe_status)==44),:)=[0 0 2].*ones(3,length(find((apoe_status)==44)))';

% read genotypes
file_path=strcat(genome_input,'.bim');
fileID = fopen(file_path, 'r');
geno_bim = textscan(fileID, '%s %s %f %d %s %s', 'Delimiter', '\t');  
fclose(fileID);
snpswithRS=geno_bim{2};
Msnp=length(snpswithRS); 

file_path=strcat(genome_input,'.fam');
fileID = fopen(file_path, 'r');
geno_fam = textscan(fileID, '%s %s %s %s %s %s');  
fclose(fileID);
geno_IID=geno_fam{2};
[C,ia,ib] = intersect(id_var,geno_IID,'legacy'); 
%assert(all(ia'==1:length(ia)))  % we assume that genotype data is available for all participants, otherwise we need to subset phenotype data accrding to IA

G1=PlinkRead_binary2(length(geno_IID),1:Msnp,genome_input);
eps2=find(strcmp(snpswithRS,'rs7412')==1);
eps4=find(strcmp(snpswithRS,'rs429358')==1);
G1(:,[min(eps2,eps4):max(eps2, eps4)])=0; % eliminate two apoe variants (we'll use APOE e2 & e4 dosages as covariates)
G=double(G1(ib, :));
clear('G1');

%%

time_var = time_var(ia,:); 
event_var = event_var(ia,:);  % 1 = cases, 0 = controls
id_var = id_var(ia,:);
country_var = country_var(ia,:);
apoe_status = apoe_status(ia,:);
apoe_status2 = apoe_status2(ia,:);


%%
options = glmnetSet;
options.dfmax=dfmax;
options.intr=false;  % silence the warning about intercept; it's not needed for coxph fit

% http://cran.nexr.com/web/packages/glmnet/vignettes/Coxnet.pdf
% y is an n length vector of failure/censoring times, and status is an n
% length vector with each entry, a 1 or a 0, indicating whether the corresponding
% entry in y is indicative of a failure time or right censoring time (1 for failure, 0
% for censoring)
y=cat(2, time_var, event_var);
[sum(event_var) sum(~event_var)]

train_idx = true(size(y, 1), 1);

dummyvar_country = dummyvar(categorical(country_var));
dummyvar_country = dummyvar_country(:, sum(dummyvar_country(train_idx, :)) ~= 0);
[~, idx] = max(sum(dummyvar_country));
dummyvar_country(:, idx) = [];  % drop most frequent country

% model specificiation: apoe and genome; control for country status
X0=[G, apoe_status2(:,[1 3])]; X=[X0, dummyvar_country]; pf=zeros(1,size(X,2)); pf(1:(size(G,2))) = 1;

if all(pf==0)
    [betas, logl, H, stats] = coxphfit(X(train_idx, :), time_var(train_idx),'censoring',1-event_var(train_idx));
else
    options.alpha = alpha;
    options.penalty_factor=pf;
    options.lambda = [];
    fit1=cvglmnet(X(train_idx, :), y(train_idx, :), 'cox', options, 'deviance', nfolds);
    [~, indlam]=min( abs(fit1.lambda - fit1.lambda_min) );
    betas=fit1.glmnet_fit.beta(:, indlam);
end

% exclude variables used during training to control for batch effects (e.g. country-specific effects) 
betas = betas(1:size(X0, 2));

%phs = X0 * betas;

tbl_betas=table([geno_bim{2}' {'apoe2'} {'apoe4'} ]', betas);
%writetable(tbl_betas,out_fname); % has header
writecell(table2cell(tbl_betas), sprintf('%s.csv', out_prefix)); % without header


