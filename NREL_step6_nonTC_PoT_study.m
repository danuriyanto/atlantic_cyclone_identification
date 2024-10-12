clc; clear; close all
addpath cyclone_id_functions/
load("alldist_all_hurdat.mat")
datasource = 'NREL';
IID_NREL = alldist_all_hurdat(5).result;
percentile = [75 80 85 90 95];
Hs50 = zeros(12,numel(percentile));
Hs500 = Hs50;
lam = Hs500;
%%
for loc = 1:12
    obsyear = 32;
    site_desc = IID_NREL(loc).name;
    extremes_timehistory = IID_NREL(loc).nonTC_data;

    %% get the nonTC distribution
    parfor th = 1:numel(percentile)
        threshold_percentile = percentile(th);
        [points_nonTC,thres] = get_PoT_data(extremes_timehistory,threshold_percentile,days(3));
        points_nonTC = sort(points_nonTC.Hs,"ascend");
        lam(loc,th) = numel(points_nonTC)/obsyear;
        % fitting to different distribution options
        dist_nonTC_GPD     = GPD_LS_or_MLE(points_nonTC,thres);
        stats_method_nonTC = stats_method(points_nonTC, obsyear, dist_nonTC_GPD,[50 500]', IID_NREL(loc).name, ['nonTC th = ' num2str(threshold_percentile)],'y');
        Hs50(loc,th)        = stats_method_nonTC.MRI_requested(1,2);
        Hs500(loc,th)       = stats_method_nonTC.MRI_requested(2,2);
        fprintf('threshold = %0.0f percent loc = %0.0f',threshold_percentile,loc)
    end

end

%% plotting
f = figure;
f.Position = [100 100 900 600];

for i = 1:12
    sitede = IID_NREL(i).name;
    subplot(4,3,i)
    plot(lam(i,:),Hs50(i,:),Marker=".",MarkerSize=10,DisplayName='Hs_{50} (m)')
    hold on
    plot(lam(i,:),Hs500(i,:),Marker=".",MarkerSize=10,DisplayName='Hs_{500} (m)')
    title(sitede)
    grid on
    ylim([0 14])
    yticks(0:2:14)
    xlim([1 10])
    xticks([1:1:10])
    if i == 12
        legend(Location="best")
    end
end
ax = axes(Visible='off');
ax.Title.Visible = 'on';
ax.XLabel.Visible = 'on';
ax.YLabel.Visible = 'on';
ax.YLabel.String = 'Hs (m)';
ax.XLabel.String = '\lambda';
ax.YLabel.FontWeight = "bold";
ax.XLabel.FontWeight = "bold";
exportgraphics(f,'PoT_nonTC_study.png',Resolution=600)