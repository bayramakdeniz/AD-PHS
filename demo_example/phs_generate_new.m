function outfile = phs_generate_new(gmat,age,casevec,varargin)


% phs_generate_new.m
%
% Purpose: Select SNPs and calculate PHS betas for gmat
% 
% Usage: phs_generate_new(gmat,age,casevec,<varargin>)
%
% Created: 03-06-17 by Tyler Seibert, Anders Dale, Chun Chieh Fan
% Copyright 03-06-17 by Creators; use only by permission
% See Change Log at end of file for modifications
% 
% Required Arguments:
%   gmat            =   [numeric] n x p array of genotype counts, where n
%                       is the number of subjects, and p is the number of
%                       SNPs (corresponds to snpidlist order). Entries
%                       should be 0,1,2 for number of allele copies (e.g.,
%                       could make it number of copies of minor allele, but
%                       this is not important, as long as all data
%                       internally consistent). Missing SNPs should have
%                       negative entries (e.g., -1) instead of 0,1,2. Call
%                       rates should already be checked, and desired QC
%                       steps should already be performed.
%   age             =   [double] vector containing age of each subject at
%                       either time of diagnosis or last follow-up in the
%                       case of controls.
%   casevec         =   [double] vector containing 1 for cases and 0 for
%                       controls (who are censored at time of observation).
%                       
% 
% Optional Arguments (varargin). Defaults in {};
%   snpidlist       =   [cell] vector containing rs IDs for SNPs in gmat.
%                       If empty, will not report SNP IDs in output.
%   covariates      =   [double] array with n rows where columns contain
%                       data of interest for covariates of user's choice.
%                       For example, 6 principal components for genotype
%                       could be included as 6 columns. Rows must match
%                       order of gmat for subjects. 
%                       This is specifically for covariates in the logistic
%                       regression step for SNP selection and is *not* 
%                       applied to the Cox regression that generates the 
%                       final PHS betas. {[]}
%   pthresh         =   [double] single value representing desired
%                       threshold for the p-values from the initial trend
%                       test. Only SNPs with p-value less than pthresh will
%                       be candidates for inclusion in the PHS model. This
%                       value is fairly arbitrary. In the PCa dataset, we
%                       used 1e-6, which yielded ~2400 candidate SNPs.
%                       Smaller data sets will need a higher pthresh.
%                       {1e-6}
%   pthresh2        =   [double] single value representing desired
%                       threshold for whether to include a SNP in the PHS
%                       model after testing the candidate SNP in a general
%                       linear model with age and <covariates>. If the
%                       p-value for adding the candidate SNP to the GLM is
%                       lower than pthresh2, the SNP is included in the Cox
%                       proportional hazards (i.e., PHS) model, and a
%                       beta-value is calculated for it. By default
%                       pthresh2 is taken to be the same as pthresh, but
%                       the user can specify a different value here. {[]}
%   outfile         =   [char] path for output file (.mat).
%                       {'new_PHS_outfile.mat'}
% 
% 
% Outputs:
%   outfile = [char] .mat file with output
% 
% Saved files:
%   outfile. .mat contents include:
%       ***
%       ***
% 
% Dependencies:
%   parse_varargin.m
%   trendTest_GWAS_amd.m
% 

defaults = {
    'snpidlist',{};
    'covariates',[];
    'pthresh',1e-6;
    'pthresh2',[];
    'outfile','new_PHS_outfile.mat';
    };

% Parse inputs:
[snpidlist,covariates,pthresh,pthresh2,outfile] = parse_varargin(varargin,defaults);


if isempty(pthresh2)
    pthresh2 = pthresh;
end
    
% build censorvec from casevec
censorvec = nan(size(casevec));
censorvec(casevec==1)=0;
censorvec(casevec==0)=1; % cases were followed to "failure." controls were censored at age of interview

% get pvec from trend test
nsnp = size(gmat,2);
pvec = NaN(1,nsnp);
for snpi = 1:nsnp
  dvec = gmat(:,snpi)>=0; % only want to use the patients who do not have missing data for this SNP
  pvec(snpi) = trendTest_GWAS_amd_Mar2017(casevec(dvec),gmat(dvec,snpi));
  if mod(snpi,1000)==0, fprintf(1,'snpi=%d of %d\n',snpi,nsnp); end
end

% sort pvec --> sv ("sorted vector")
[sv,si] = sort(pvec,'ascend');

snplist = []; % SNPs to include in final model; initialize here
num_candidate_SNPs = find(sv<pthresh,1,'last'); % faster than 1:max(find(sv<pthresh))
for topi = 1:num_candidate_SNPs
  X = double(gmat(:,[si(topi) snplist])); % add current test snp to gmat
  
  % impute missing genotype data (replace with mean):
  M = zeros(size(X));
  for mm=1:size(X,2); 
      mvec=X(:,mm); 
      X(mvec<0,mm)=mean(mvec(mvec>=0)); 
      M(mvec<0,mm)=-1;
  end
  
  % GLM to select SNPs for PHS model
  mdl = fitglm([X age covariates],casevec); % important to have X first or "2" in next line will be wrong
  pval = mdl.Coefficients.pValue(2); % model p-value for 
  if pval<pthresh2 % if under threshold, then include this SNP in PHS
    snplist = cat(2,si(topi),snplist); % add candidate SNP to list for PHS
    svec = mdl.predict;

    defvec = isfinite(age); 
%    [b logl H stats_cox] = coxphfit(X(defvec,:),age(defvec),'censoring',censorvec(defvec),'baseline');
    [b logl H stats_cox] = coxphfit(X(defvec,:),age(defvec),'censoring',censorvec(defvec)); % May want to subtract mean of each column on X for general population, and use 'baseline' option
    Hcox = H; bcox = b; % save b and H from coxphfit
    logH_pred = X*b; 
    fprintf(1, 'added SNP number %d of %d candidates\n',topi,num_candidate_SNPs);
  end
end

% clean up:
if isempty(topi)
    error('%s: final model did not include ANY SNPs\n',mfilename);
elseif topi==num_candidate_SNPs && pval<pthresh2
    fprintf('%s: final model included ALL candidate SNPs identified in trend test\n',mfilename);
elseif pval>=pthresh2 % if last candidate SNP tested was not included in final model (this is the norm)
    X = X(:,2:end); % remove last candidate SNP
end
if numel(X(1,:))~=numel(snplist)
    error('%s: Something is wrong. Number of SNPs in X and snplist should match\n',mfilename);
end

modelsnps = snpidlist(snplist); % names of final model snps

% Save output
clear gmat mdl I % too big
save(outfile,'-v7.3')
fprintf('%s: wrote %s\n',mfilename,outfile);


% ***** To Do: 
% -create PHS_plots and/or PHS_figures for here
%   --consider re-adding the AUCs...
%   --get the colororder and everything back
% -define which variables really need to be saved here
% -finish documentation above
% 


