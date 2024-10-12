clc; clear; close all
addpath cyclone_id_functions/
dir_stats_method = dir("stats_method*mat");
dir_alldist      = dir("alldist*mat");
method_status    = {'all hurdat' 'hurricane only' 'no ETC'};
sites = readtable("data/12_sites_final.csv");
for loc=1:12
    for i=1:numel(dir_stats_method)
        status_filter       = dir_stats_method(i).name;
        status_filter       = split(status_filter,'_');
        load(dir_stats_method(i).name)
        load(dir_alldist(i).name)
        IID_NREL = alldist_all_hurdat(5).result;

        % load the fitted NREL 2-params Weibull
        MRP_fitted_TC       = stats_method_TC(loc,4).MRP_fitted;
        MRP_fitted_nonTC    = stats_method_nonTC(loc,1).MRP_fitted;

        % fitted distribution
        dist_TC             = stats_method_TC(loc,4).dist_object;
        dist_nonTC          = stats_method_nonTC(loc,1).dist_object;

        % calculate the TC and nonTC arrival rate
        lambda_TC           = height(IID_NREL(loc).MIS_TC)/32;
        lambda_nonTC        = 1;
        RP                  = [50 500]';
        empiricalResult     = empiricalMRP(lambda_TC,lambda_nonTC,dist_TC,dist_nonTC,RP);
        combined_table      = empiricalResult.MRI_table;

        % MRP
        MRI_requested_TC{loc,i}            = stats_method_TC(loc,4).MRI_requested;
        MRI_requested_nonTC{loc,i}         = stats_method_nonTC(loc,1).MRI_requested;
        MRI_requested_combined{loc,i}      = empiricalResult.RP_request{:,:};


        % combine without separation
        extremes_timehistory = IID_NREL(loc).Hs_all;
        [extremes_TT,extreme_vectors,dist_combined_ns] = ...
        rlargest_CDF(extremes_timehistory, 1, sites.Location{loc},"n","Weibull",[50 500]);
        stats_method_combined_ns = stats_method(extreme_vectors, 32, dist_combined_ns, [50 500]', sites.Location{loc}, 'combined no-separation',"n");
        stats_method_combined_ns_all(loc,i) = stats_method_combined_ns;

        %% plot MRP
        strt1 = find(MRP_fitted_nonTC(:,1)==1.3);
        ff    = find(combined_table.T>1.3);
        strt2 = ff(1);
        f = figure(Visible="off");
        f.Position = [476 360 560 300];
        semilogx(MRP_fitted_TC(strt1:end,1),MRP_fitted_TC(strt1:end,2),      DisplayName='TC fitting - 2p Weibull', ...
            Color=[0.4660 0.6740 0.1880],LineWidth=1.5)
        hold on
        semilogx(MRP_fitted_nonTC(strt1:end,1),MRP_fitted_nonTC(strt1:end,2),DisplayName='nonTC fitting - GEV', ...
           Color='r',LineWidth=1.5)
        semilogx(combined_table.T(strt2:end),combined_table.MRI(strt2:end),  DisplayName='combined TC & nonTC distribution', ...
            Color='b',LineWidth=1.5)
        semilogx(stats_method_combined_ns.MRP_fitted(strt1:end,1), ...
            stats_method_combined_ns.MRP_fitted(strt1:end,2),                DisplayName='no separation - 2p Weibull',...
            Color='m',LineWidth=1.5)
        box on
        grid on
        xlim([1 10000])
        fontsize(f,14,"points")
        legend(Location="northwest",FontSize=10)
        xlabel('MRP (years)',FontSize=14)
        ylabel('Hs (m)',FontSize=14)
        % title([IID_NREL(loc).name ' ' method_status{i}])
        exportgraphics(f,['weibull_2p_filtering_comparison/' method_status{i} '_' IID_NREL(loc).name '.png'],Resolution=450)
    end
end

sitename = sites.Location;
TC = [vertcat(MRI_requested_TC{:,2}) vertcat(MRI_requested_TC{:,1}) vertcat(MRI_requested_TC{:,3})];
TC(:,[3 5]) = [];
TC = array2table(TC,VariableNames={'RP' 'only TD TS HU' 'excluding EX'  'all HURDAT'});
TC50 = TC(TC.RP==50,:);
TC50 = addvars(TC50,sitename,Before='RP');
TC500 = TC(TC.RP==500,:);
TC500 = addvars(TC500,sitename,Before='RP');


nonTC = [vertcat(MRI_requested_nonTC{:,2}) vertcat(MRI_requested_nonTC{:,1}) vertcat(MRI_requested_nonTC{:,3})];
nonTC(:,[3 5]) = [];
nonTC = array2table(nonTC,VariableNames={'RP' 'only TD TS HU' 'excluding EX'  'all HURDAT'});
nonTC50 = nonTC(nonTC.RP==50,:);
nonTC50 = addvars(nonTC50,sitename,Before='RP');
nonTC500 = nonTC(nonTC.RP==500,:);
nonTC500 = addvars(nonTC500,sitename,Before='RP');

combined = [vertcat(MRI_requested_combined{:,2}) vertcat(MRI_requested_combined{:,1}) vertcat(MRI_requested_combined{:,3})];
combined(:,[3 5]) = [];
combined = array2table(combined,VariableNames={'RP' 'only TD TS HU' 'excluding EX'  'all HURDAT'});
combined50 = combined(combined.RP==50,:);
combined50 = addvars(combined50,sitename,Before='RP');
combined500 = combined(combined.RP==500,:);
combined500 = addvars(combined500,sitename,Before='RP');

combined_ns = [vertcat(stats_method_combined_ns_all(:,1).MRI_requested) ...
               vertcat(stats_method_combined_ns_all(:,2).MRI_requested) ...
               vertcat(stats_method_combined_ns_all(:,3).MRI_requested)];
combined_ns(:,[3 5]) = [];
combined_ns = array2table(combined_ns,VariableNames={'RP' 'only TD TS HU' 'excluding EX'  'all HURDAT'});
combined_ns50 = combined_ns(combined_ns.RP==50,:);
combined_ns50 = addvars(combined_ns50,sitename,Before='RP');
combined_ns500 = combined_ns(combined_ns.RP==500,:);
combined_ns500 = addvars(combined_ns500,sitename,Before='RP');

all_50 = table(sitename,TC50.("all HURDAT"),nonTC50.("all HURDAT"),combined50.("all HURDAT"),combined_ns50.("all HURDAT"), ...
    VariableNames={'sitename' 'TC' 'nonTC' 'combined' 'no separation'})
all_500 = table(sitename,TC500.("all HURDAT"),nonTC500.("all HURDAT"),combined500.("all HURDAT"),combined_ns500.("all HURDAT"), ...
    VariableNames={'sitename' 'TC' 'nonTC' 'combined' 'no separation'})


writetable(TC50,'MRI_result_table/TC50_raw.csv')
writetable(TC500,'MRI_result_table/TC500_raw.csv')

writetable(nonTC50,'MRI_result_table/nonTC50_raw.csv')
writetable(nonTC500,'MRI_result_table/nonTC500_raw.csv')

writetable(combined50,'MRI_result_table/combined50_raw.csv')
writetable(combined500,'MRI_result_table/combined500_raw.csv')

writetable(combined_ns50,'MRI_result_table/combined50ns_raw.csv')
writetable(combined_ns500,'MRI_result_table/combined500ns_raw.csv')
