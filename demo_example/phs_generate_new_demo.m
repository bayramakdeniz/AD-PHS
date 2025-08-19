

% generate new PHS model from genotype / age data:
load new_PHS_input_example % some junk data I made for testing
pthresh = 1e-2; 
outfile = phs_generate_new(q.gmat,q.age,q.casevec,'pthresh',pthresh,'snpidlist',q.snpidlist)


% define population absolute incidence and baseline hazard functions:
agevals = 40:100;
a = 0.0700/100; b = 0.0753; pcainc = a*exp(b*(agevals-40)); % this is the population incidence rate for PCa, as described in our paper
H0 = a/b*exp(b*(agevals-40)); % compute baseline hazard
H0 = H0-H0(1); % assume cumulative hazard == 0 before age 40


% plot figures for above model:
load(outfile);
PHSvec = logH_pred; % just renaming for clarity
PHS_ref_population = logH_pred(age<70 & casevec==0); % controls <70
savefig = 0;
groupname = 'newPHSinputexample'; % just a name for figure titles to indicate dataset
diseasename = 'PCa'; % just a name fof figure titles
N.discoveryset = phs_plots(PHSvec,age,censorvec,agevals,H0,'reference_PHS',PHS_ref_population,'groupname',groupname,'diseasename',diseasename,'savefig',savefig); % Plot figures



