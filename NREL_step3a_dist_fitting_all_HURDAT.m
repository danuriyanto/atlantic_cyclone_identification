clc; clear; close all
addpath cyclone_id_functions/
load("alldist_all_hurdat.mat")
datasource = 'NREL';
IID_NREL = alldist_all_hurdat(5).result;
for loc = 1:12
    clearvars maxHs
    obsyear = 32;
    site_desc = IID_NREL(loc).name;

    %% TC fitting
    maxHs = IID_NREL(loc).MIS_TC;
    points_TC = sort(maxHs.Hs,"ascend");

    %% fitting to distributions
    dist_TC_GEV      = GEV_LS_or_MLE(points_TC);
    dist_TC_gumbel   = fitdist(points_TC,"ExtremeValue");
    dist_TC_logn     = fitdist(points_TC,"Lognormal");
    dist_TC_weibull  = fitdist(points_TC,"Weibull");

    % create a struct to curate all the fittings
    dist_TC = struct();
    dist_TC(1).name = 'GEV with bounded k';
    dist_TC(1).dist = dist_TC_GEV;
    dist_TC(2).name = 'Gumbel';
    dist_TC(2).dist = dist_TC_gumbel;
    dist_TC(3).name = 'Lognormal';
    dist_TC(3).dist = dist_TC_logn;
    dist_TC(4).name = '2-params Weibull';
    dist_TC(4).dist = dist_TC_weibull;

    %% get the nonTC distribution

    % get the IID using yearly maxima (r = 1)
    points_nonTC = sort(IID_NREL(loc).nonTC_AM.Hs,"ascend");
    % fitting to different distribution options
    dist_nonTC_GEV      = GEV_LS_or_MLE(points_nonTC);
    dist_nonTC_gumbel   = fitdist(points_nonTC,"ExtremeValue");
    dist_nonTC_logn     = fitdist(points_nonTC,"Lognormal");
    dist_nonTC_weibull  = fitdist(points_nonTC,"Weibull");

    % create a struct to curate all the fittings
    dist_nonTC = struct();
    dist_nonTC(1).name   = 'GEV with bounded k';
    dist_nonTC(1).dist   = dist_nonTC_GEV;
    dist_nonTC(2).name   = 'Gumbel';
    dist_nonTC(2).dist   = dist_nonTC_gumbel;
    dist_nonTC(3).name   = 'Lognormal';
    dist_nonTC(3).dist   = dist_nonTC_logn;
    dist_nonTC(4).name   = '2-params Weibull';
    dist_nonTC(4).dist   = dist_nonTC_weibull;

    %% set 2: statistical method to calculate the MRI
    % calculate the diagnostics
    for i=1:numel(dist_TC)
        stats_method_TC(loc,i)         = stats_method(points_TC, obsyear, dist_TC(i).dist,[50 500]', IID_NREL(loc).name, ['TC ' dist_TC(i).name],'y');
        stats_method_nonTC(loc,i)      = stats_method(points_nonTC, obsyear, dist_nonTC(i).dist,[50 500]', IID_NREL(loc).name, ['nonTC ' dist_nonTC(i).name],'y');
    end
    close all

    save('stats_method_fitting_all_HURDAT',"stats_method_TC","stats_method_nonTC")
    
end

%% compare all stats method with the Surrogate
close all
load("../DATASET_SURROGATE_NEW/MRI_hurricane_surrogate_peaks.mat")

for s=1:12

    %% plot the TC comparison
    f_TC = figure;
    f_TC.Position = [100 100 1120 420];

    cc = {
    [0, 0, 0],           % Black
    [139, 0, 0]/255,     % Dark Red
    [128, 128, 0]/255,   % Olive (Dark Yellow)
    [0, 0, 139]/255,     % Dark Blue
    [255, 140, 0]/255,   % Dark Orange
    [128, 0, 128]/255    % Purple
        };

        
    subplot(121)
    semilogx(MRP_hurricane_surrogate_peaks(s).RP,MRP_hurricane_surrogate_peaks(s).MRI,...
        Marker=".",LineStyle="-",MarkerEdgeColor=cc{1},Color=cc{1},...
        DisplayName='100,000yr catalog')
    hold on

    semilogx(stats_method_TC(s,2).MRP_empirical(:,1),stats_method_TC(s,2).MRP_empirical(:,2),...
        Marker="o",LineStyle="--",MarkerEdgeColor=cc{1},Color=cc{1},LineWidth=0.7,...
        DisplayName='empirical MRP')


    for i=1:4
        dispnames = split(stats_method_TC(s,i).dist_type,'TC ');
        if stats_method_TC(s,i).ad_test == 0 && stats_method_TC(s,i).ks_test == 0
        semilogx(stats_method_TC(s,i).MRP_fitted(:,1),stats_method_TC(s,i).MRP_fitted(:,2),...
             LineStyle="-",MarkerEdgeColor=cc{i+2},Color=cc{i+2},LineWidth=1.2,...
             DisplayName=dispnames{2})
        end
    end
    grid on
    fontsize(f_TC,14,"points")
    legend(Location="southeast",FontSize=14)
    xlabel('MRP (years)',FontSize=16)
    ylabel('Hs (m)',FontSize=16)
    title([IID_NREL(s).name ' TC (Method of Independent Storm)'])
    ylim([0 17])
    xlim([1 1e4])
    

    %% plot the nonTC comparison
    % f_nonTC = figure;
    % f_nonTC.Position = [100 100 1120/2 420];
    sst = 161;
    subplot(122)
    semilogx(stats_method_nonTC(s,2).MRP_empirical(:,1),stats_method_nonTC(s,2).MRP_empirical(:,2),...
        Marker="o",LineStyle="--",MarkerEdgeColor=cc{1},Color=cc{1},LineWidth=0.7,...
        DisplayName='empirical MRP')
    hold on

    for i=1:4
        dispnames = split(stats_method_nonTC(s,i).dist_type,'nonTC ');
        if stats_method_nonTC(s,i).ad_test == 0 && stats_method_nonTC(s,i).ks_test == 0
        semilogx(stats_method_nonTC(s,i).MRP_fitted(sst:end,1),stats_method_nonTC(s,i).MRP_fitted(sst:end,2),...
             LineStyle="-",MarkerEdgeColor=cc{i+2},Color=cc{i+2},LineWidth=1.2,...
             DisplayName=dispnames{2})
        end
    end
    grid on
    fontsize(f_TC,14,"points")
    legend(Location="southeast",FontSize=14)
    xlabel('MRP (years)',FontSize=16)
    ylabel('Hs (m)',FontSize=16)
    title([IID_NREL(s).name ' nonTC (Annual Maxima)'])
    ylim([0 17])
    xlim([1 1e4])
    
     exportgraphics(f_TC,['method_comparison_plots/'...
        'all_HURDAT_' IID_NREL(s).name ...
        '_method_comp.png'],Resolution=300)
end