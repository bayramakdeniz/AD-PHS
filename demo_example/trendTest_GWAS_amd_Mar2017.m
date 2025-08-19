function [pvec zvec] = trendTest_GWAS_amd(yvec,gmat)
%function [pvec zvec] = trendTest_GWAS_amd(yvec,gmat)

% Created by Anders Dale. Copied current version to here in March 2017.

nsnp  = size(gmat,2);

zvec    = NaN(nsnp,1);
logpvec = NaN(nsnp,1);
pvec    = NaN(nsnp,1);

for i=1:nsnp
  s0 = sum(gmat(yvec==0,i)==0);
  s1 = sum(gmat(yvec==0,i)==1);
  s2 = sum(gmat(yvec==0,i)==2);
  nsubc = s0+s1+s2;

  r0 = sum(gmat(yvec==1,i)==0);
  r1 = sum(gmat(yvec==1,i)==1);
  r2 = sum(gmat(yvec==1,i)==2);
  nsuba = r0+r1+r2;

  nsub = nsubc+nsuba;
  
  n1 = s1+r1;
  n2 = s2+r2;
  
  xBarC = (2*s2 + s1)/nsubc;
  xBarA = (2*r2 + r1)/nsuba;
  xBar  = (2*n2 + n1)/nsub;
  
  z = (xBarC - xBarA)/sqrt( (4*n2 + n1 - nsub*(xBar^2))/(nsubc * nsuba) );
  p = 2*normcdf(-abs(z));
  
  zvec(i)    = z;  
  pvec(i)    = p;
%  logpvec(i) = -log10(p);
end

%keyboard

% Fixed:
%   Remove assumption of sorted cases and controls
%   Handle missing genotypes

% ToDo
%   Reverse sign of zvec?
