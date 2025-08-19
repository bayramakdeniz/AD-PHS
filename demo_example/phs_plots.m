function N = phs_plots(PHSvec,age,censorvec,agevals,H0,varargin)

% phs_plots.m
%
% Purpose: Standard plots for PHS data
% 
% Usage: phs_plots(PHSvec,age,censorvec,agevals,H0,<varargin>)
%
% Created: 03-16-17 by Tyler Seibert, Anders Dale
% Copyright 03-16-17 by Creators; use only by permission
% See Change Log at end of file for modifications
% 
% Required Arguments:
%   PHSvec          =   [double] n x 1 array of calculated PHS for each of
%                       n subjects in dataset.
%   age             =   [double] vector containing age of each subject at
%                       either time of diagnosis or last follow-up for 
%                       controls.
%   censorvec       =   [numeric] n x 1 vector containing 1 for each
%                       subject to be censored and 0 for others (typically
%                       1 for controls and 0 for cases). 
%   agevals         =   [double] age vector of interest. Typically want to
%                       have it go a bit older than you want to plot
%                       because of the odd effects at high range of plot. 
%                       For example, we use 100 as max and then only plot 
%                       to 95. {[40:100]}
%   H0              =   [double] vector of baseline hazard for all values
%                       in agevals. Length must match agevals. 
% 
% Optional Arguments (varargin). Defaults in {};
%   reference_PHS   =   [double] vector of calculated PHS for some
%                       population of subjects that will be defined as the
%                       reference for percentiles of PHS. E.g., the
%                       99th%ile for PHS will be set at the 99th%ile of
%                       this population. If not defined, will just use
%                       PHSvec.
%   groupname       =   [char] name of dataset for plot titles
%                       {'newPHSdata'}
%   fignum          =   [double] Starting index for figures. {0}
%   plotideal       =   [0|1] toggle plotting K-M estimates {1}
%   plotempirical   =   [0|1] toggle plotting empirical survival data {1}
%   savefig         =   [0|1] toggle saving figures {0}
%   prctile_list    =   [double] list of percentiles desired for plots. By
%                       default, will use [1 5 20 50 80 95 99]; Note that
%                       if this is changed, colororder should also be
%                       appropriately modified. 
%   colororder      =   [double] by default, calls phs_colorder.m
%   prctile_ranges  =   [double] ranges for plots for empirical plots. If
%                       default not used, will need to modify kmcolororder
%                       appropriately, as well. 
%   kmcolororder    =   [double] by default, takes subset of colororder
%   diseasename     =   [char] name of condition studied for plot titles
%   prefix          =   [char] prefix for saved figure names
%                       {'newPHS_plot_fig'}
% 
% 
% Outputs:
%   N               =   [table] table of ***
% 
% Saved files:
%   K-M plots as requested with <plotideal> and <plotempirical>, with names
%   according to <prefix>.
% 
% Dependencies:
%   parse_varargin.m
%   phs_colororder.m
% 


defaults = {
    'reference_PHS',[]; % use all PHSvec
    'groupname','newPHSdata';
    'fignum',0;
    'plotideal',1;
    'plotempirical',1;
    'savefig',0;
    'prctile_list',[1 5 20 50 80 95 99];
    'colororder',[]; % get from phs_colorder. If numel prctile_list customized, will need to adjust
    'prctile_ranges',[0 20; 30 70; 80 98; 98 100]; % for K-M plots: "<20", "50", ">80", ">98"; We have also done: [0 2; 2 10; 10 30; 30 70; 70 90; 90 98; 98 100];
    'kmcolororder',[]; % colors for prctile_ranges; by default
    'diseasename','Disease';
    'prefix','newPHS_plot_fig';
    };

    
% Parse inputs:
[reference_PHS,groupname,fignum,plotideal,plotempirical,savefig,prctile_list,colororder,prctile_ranges,kmcolororder,diseasename,prefix] = parse_varargin(varargin,defaults);

if isempty(reference_PHS)
    reference_PHS = PHSvec; 
    fprintf('%s: no reference (population) PHS vector provided so using all values in PHSvec as reference for percentiles\n',mfilename);
end
if isempty(colororder)
    colororder = phs_colororder;
end

% get linestyle and color lists:
linestylelist = cell(1,numel(prctile_list));
linestylelist(:) = {'-'};
linestylelist(prctile_list==50) = {'--'};


% define percentiles and x vector:
critvalvec = prctile(reference_PHS,prctile_list); 
x = agevals; % just trying to preserve the old code for comparisons

% plot ideal K-M estimates:
if plotideal
    figure(fignum+11); clf; hold on; 
    title([diseasename '-free survival, by PHS percentile']); 
    h=xlabel('Age'); set(h,'FontSize',12,'FontWeight','bold'); h=ylabel('Proportion not diagnosed with prostate cancer'); set(h,'FontSize',12,'FontWeight','bold');

    figure(fignum+12); clf; hold on; 
    title(['Incidence rate, by PHS percentile']);
    h=xlabel('Age'); set(h,'FontSize',12,'FontWeight','bold'); h=ylabel('Incidence rate (%)'); set(h,'FontSize',12,'FontWeight','bold');

    for cnt = 1:length(critvalvec)
      H = H0*exp(critvalvec(cnt)-critvalvec(prctile_list==50)); S = exp(-H); % calculate hazard for this PHS percentile
      
      figure(fignum+11); plot(x,S,'LineWidth',3,'color',colororder(cnt,:),'LineStyle',linestylelist{cnt}); 
      xlim([min(agevals) max(agevals)-5]); ylim([0 1]);
      
      figure(fignum+12); plot(x,gradient(H,x(2)-x(1))*100,'LineWidth',3,'color',colororder(cnt,:),'LineStyle',linestylelist{cnt}); 
      xlim([min(agevals) max(agevals)-5]);
      
      legends{cnt} = sprintf('PHS percentile %2d',prctile_list(cnt));
    end

    fignum1=fignum+11; figure(fignum1); h=legend(legends,'Location','SW'); set(h,'FontSize',12,'FontWeight','bold');
    if savefig; print([prefix num2str(fignum1)],'-djpeg99','-r500'); end % save 501dpi jpeg with 99% quality

    fignum1=fignum+12; figure(fignum1); h=legend(legends,'Location','NW'); set(h,'FontSize',12,'FontWeight','bold');
    if savefig; print([prefix num2str(fignum1)],'-djpeg99','-r500'); end % save 501dpi jpeg with 99% quality

else
    for cnt = 1:length(critvalvec)
        legends{cnt} = sprintf('PHS percentile %2d',prctile_list(cnt));
    end
end
    



% Kaplan-Meier plots for actual data (empirical):
if plotempirical
    
    fignum1=fignum+14; figure(fignum1); figure(fignum1);
    [f1,x1,f1_low,f1_up] = ecdf(age,'censoring',censorvec,'function','survivor');
    ax1=stairs(x1,f1); ax1.LineWidth = 3; hold on; stairs(x1,f1_low,':'); stairs(x1,f1_up,':');
    ylim([0 1]);
    title(['Kaplan-Meier curve for all ' groupname ' patients']);
    if savefig; print([prefix num2str(fignum1)],'-djpeg99','-r500'); end % save 501dpi jpeg with 99% quality

    % K-M plots for percentile ranges
    titlestr = sprintf('Kaplan-Meier estimates of %s-free survival for %s, by PHS percentile',diseasename,groupname);
    fignum1=fignum+15; figure(fignum1); hold on; h44=title(titlestr); 
    set(h44,'FontSize',12,'FontWeight','bold');
    h=xlabel('Age'); set(h,'FontSize',12,'FontWeight','bold'); h=ylabel('Probability of PCa-free survival'); set(h,'FontSize',12,'FontWeight','bold');
    legends_handles=[];tmplegends={};
    if isempty(kmcolororder)
        kmcolororder = colororder([3:5 7],:);
        kmcolororder(2,:) = [40 40 42]./255;
    end
    for rr = 1:numel(prctile_ranges(:,1));
        tmprr = prctile(reference_PHS,prctile_ranges(rr,:));
        inds = find(PHSvec>tmprr(1) & PHSvec<tmprr(2));
        N(rr) = numel(inds);
        [f1,x1,f1_low,f1_up] = ecdf(age(inds),'censoring',censorvec(inds),'function','survivor');

        figure(fignum1); ax1 = stairs(x1,f1); %hold on; stairs(x1,f1_low,':'); stairs(x1,f1_up,':');
        ax1.Color = kmcolororder(rr,:); ax1.LineWidth = 3;
        ax2=stairs(x1,f1_low,':'); ax3=stairs(x1,f1_up,':'); ax2.Color=ax1.Color; ax3.Color=ax1.Color; ax2.LineWidth=1; ax3.LineWidth=1;
        tmplegends{rr} = sprintf('PHS percentile %d-%d',prctile_ranges(rr,:));
        legends_handles(rr) = ax1;
        tablename=sprintf('%s_%s_prct%d-%d.csv',mfilename,groupname,prctile_ranges(rr,1),prctile_ranges(rr,2));
        z=table(age(inds),censorvec(inds),'VariableNames',{'age','censorvec'});
        if savefig; writetable(z,tablename); end
    end
    prctile_ranges
    h1=legend(legends_handles,tmplegends,'Location','SW'); set(h1,'FontSize',12,'FontWeight','bold'); clear tmplegends
    fprintf('%s: Note error in legend due to confidence intervals. Modify in Photoshop, etc.\n',mfilename);
    xlim([min(agevals) max(agevals)-5]);
    ylim([0 1]);

    if savefig; print([prefix num2str(fignum1) '_withCIs'],'-djpeg99','-r500'); end % save 501dpi jpeg with 99% quality
    % note: I put a breakpoint on the savefig line and manually moved the title upward before the save
end

% CHANGE LOG:
% Created 03-16-17 by Tyler Seibert
% 03-27-17 --> changed ylim for the empirical survival plots to [0 1]
% 
% 
