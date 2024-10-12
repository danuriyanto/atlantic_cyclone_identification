clc; clear; close all
addpath cyclone_id_functions/
load("alldist_all_hurdat.mat")
datasource = 'NREL';
IID_NREL = alldist_all_hurdat(5).result;
Hs50 = zeros(12,10);
Hs500 = Hs50;
%%
for loc = 1:12
    obsyear = 32;
    site_desc = IID_NREL(loc).name;
    extremes_timehistory = IID_NREL(loc).nonTC_data;

    %% get the nonTC distribution
    parfor r = 1:10
        [~,points_nonTC,~] = rlargest_CDF(extremes_timehistory, r, site_desc,"n","GEV",[50 500]);
        points_nonTC = sort(points_nonTC,"ascend");
        % fitting to different distribution options
        dist_nonTC_GEV     = GEV_LS_or_MLE(points_nonTC);
        stats_method_nonTC = stats_method(points_nonTC, obsyear, dist_nonTC_GEV,[50 500]', IID_NREL(loc).name, ['nonTC r = ' num2str(r)],'y');
        Hs50(loc,r)        = stats_method_nonTC.MRI_requested(1,2);
        Hs500(loc,r)       = stats_method_nonTC.MRI_requested(2,2);
        fprintf('r = %0.0f loc = %0.0f',r,loc)
    end
    
end

%% plotting
f = figure;
f.Position = [100 100 900 600];
r_array = 1:1:10;
for i = 1:12
    sitede = IID_NREL(i).name;
    subplot(4,3,i)
    plot(r_array,Hs50(i,:),Marker=".",MarkerSize=10,DisplayName='Hs_{50} (m)')
    hold on
    plot(r_array,Hs500(i,:),Marker=".",MarkerSize=10,DisplayName='Hs_{500} (m)')
    title(sitede)
    grid on
    ylim([0 14])
    yticks(0:2:14)
    xlim([1 10])
    xticks([0:1:10])
    if i == 12
        legend(Location="northeast")
    end
end
ax = axes(Visible='off');
ax.Title.Visible = 'on';
ax.XLabel.Visible = 'on';
ax.YLabel.Visible = 'on';
ax.YLabel.String = 'Hs (m)';
ax.XLabel.String = 'R - yearly Maxima (\lambda)';
ax.YLabel.FontWeight = "bold";
ax.XLabel.FontWeight = "bold";
exportgraphics(f,'r_yearly_max_nonTC_study.png',Resolution=600)