% This is the matlab code section we used for plos medicine paper

X2 = gmat_test;
X2_pop_cen_b = X2_pop*bmean - pop_mean2'*bmean;
X2_cen_b = X2*bmean - pop_mean2'*bmean;
Xq = quantile(X2_cen_b(isfinite(X2_cen_b)), [0.01 0.1 0.30 0.4 0.6 0.7 0.9 0.99]);
Xq = [-Inf Xq Inf];
for i = 1:5;
  indimat(:,i) = (X2_cen_b > Xq(1+2*(i-1))).*(X2_cen_b <= Xq(2*i));
end
indimat = logical(indimat);
[b,logl,H,stats] = coxphfit(X2_cen_b, T,'censoring',censorvec,'baseline',0);
xvec = H(:,1); fvec = H(:,2);
[xx logf0] = smoothfun1d(xvec,log(fvec),agevals,1e5); f0 = exp(logf0);
Xq2 = quantile(X2_cen_b(isfinite(X2_cen_b)), [0.25 0.5 0.75]);
Xq2 = [-Inf Xq2 Inf];
for i = 1:4;
  indimat2(:,i) = (X2_cen_b > Xq2(i)).*(X2_cen_b <= Xq2(i+1));
end
indimat2(:,5) = (X2_cen_b > Xq2(2)).*(X2_cen_b <= Xq2(4));
indimat2 = logical(indimat2);
figure(8);clf;hold on;
cseq = [1,5,3,2,8];
ColorOrder = get(gca,'ColorOrder');
ColorOrder(8,:) = [0.5 0.5 0.5];
for i = 1:5;
  [f,x,flow,fup] = ecdf(T(indimat2(:,i)), 'censoring', censorvec(indimat2(:,i)), 'function', 'survivor');
  ax1(i) = plot(x,f, 'LineWidth',1,'color',ColorOrder(cseq(i),:)*0.6,'LineStyle', ':');
  plotshaded(x(2:end-1)', cat(2, flow(2:end-1), fup(2:end-1))',ColorOrder(cseq(i),:))
  %stairs(x,flow, 'LineWidth',1.5,'color',ColorOrder(cseq(i),:)*0.6,'LineStyle', ':');
  %stairs(x,fup, 'LineWidth',1.5,'color',ColorOrder(cseq(i),:)*0.6,'LineStyle', ':');
  M = mean(X2_cen_b(indimat2(:,i)));
  Ss = exp(-(f0*exp(M)));
  ax2(i) = plot(xx, Ss, 'LineWidth',3,'color',ColorOrder(cseq(i),:));
end
xlabel('Age');
ylabel('Survival');
axis([60 95 0 1])
legend([ax2(1) ax2(2) ax2(3) ax2(4) ax2(5)], {'1st PHS Quartile', '2nd PHS Quartile', '3rd PHS Quartile', '4th PHS Quartile', 'Uncorrected baseline'}, 'Location', 'SW');
xlabel('Age');
ylabel('Survival Proportion in ADGC phase 1');
axis([60 100 0 1])
print('-dpng','-r300','/space/syn03/1/data/cfan/ADGC/tmp/Surv_IGAP_v2_boot_ADGC_phase1_graph_compare_Cox_KM.png');
