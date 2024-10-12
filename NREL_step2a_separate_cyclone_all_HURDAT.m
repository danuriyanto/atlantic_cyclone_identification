clc; clear; close all;
addpath cyclone_id_functions/
dir_cyclone = dir("TC_*.mat");
load("data/Hs_NREL.mat")
load("data/V10_NCEP.mat")
sites = readtable("data/12_sites_final.csv");

IID_NREL_all_HURDAT = struct();
dist_array = [100 200 300 400 500 600 700 800 900 1000];
for dista = 1:numel(dist_array)
    distance_filter = dist_array(dista);

    for loc = 1:12
        clearvars -except alldist_all_hurdat dist_array dista IID_NREL_all_HURDAT loc dir_cyclone sites loc Hs_NREL distance_filter Vdata
        load(dir_cyclone(loc).name)

        % load Hs
        Hs_data = Hs_NREL(:,loc);
        Hs_data.Properties.VariableNames = {'Hs'};
        id_neg = Hs_data.Hs<=0;
        Hs_data(id_neg,:) = [];

        % load u10
        u10_all = Vdata(:,loc);
        u10_all.Properties.VariableNames = {'u10_local'};

        site_lat = sites.Latitude(loc);
        site_lon = sites.Longitude(loc)+360;
        site_desc = sites.Location{loc};
        plot_lim = [site_lat-5 site_lat+5
            site_lon-5 site_lon+5];

        %% get TC peaks
        % Initialize variables to store the results of the retiming
        TC_peaks = timetable();
        unq_TC_ID = unique(TC_data.storm_ID);
        for i = 1:numel(unq_TC_ID)
            TC_at_i = TC_data(TC_data.storm_ID == unq_TC_ID(i),:);
            [~,c0] = max(TC_at_i.Hs);
            peakstorm = TC_at_i(c0,:);
            TC_peaks = [TC_peaks; peakstorm];
        end

        % findEX = find(TC_peaks.status == "EX");
        % TC_peaks(findEX,:) = [];
        % clear out the Hs peaks where the storm peak is outside search radius
        longdist_TC             = find(TC_peaks.distance_from_position > distance_filter);
        TC_peaks(longdist_TC,:) = [];
        filter_TC               = ismember(TC_data.storm_ID,TC_peaks.storm_ID);
        TC_data_filtered        = TC_data(filter_TC,:);


        %% get nonTC extremes
        idnontc = ismember(Hs_data.Time,TC_data_filtered.Time);
        nonTC_data = Hs_data(~idnontc,:);
        r = 1;
        [extremes_TT,extreme_vectors,distparameter] = rlargest_CDF(nonTC_data, r, site_desc,"n","Weibull",[50 500]);


        %% plot the timeseries
        fig_timehistory = figure;
        fig_timehistory.Position = [100 100 800 300];
        s1 = plot(Hs_data.Time,Hs_data.Hs,DisplayName='3-hourly data', Color=[0.5 0.5 0.5], LineWidth=0.4);
        hold on

        % Define the range of years in your data
        years = year(Hs_data.Time(1)):year(Hs_data.Time(end));

        % Add subtle gray dashed lines between the years
        for y = years
            xline(datetime(y, 1, 1), '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
        end
        s2 = scatter(TC_data_filtered.Time,TC_data_filtered.Hs,Marker=".",DisplayName="TC",SizeData=20,MarkerEdgeColor="blue",LineWidth=1.1);
        s3 = scatter(TC_peaks.Time,TC_peaks.Hs,Marker="x",DisplayName="TC peaks",SizeData=20,MarkerEdgeColor="blue",LineWidth=1.5);
        s4 = scatter(extremes_TT.Time,extremes_TT.Hs,Marker="o",DisplayName="nonTC peaks",SizeData=20,MarkerEdgeColor="red",LineWidth=1.1);
        legend([s1 s3 s4], Location="best")
        % title([sites.Location{loc} ' Hs time history R = ' num2str(distance_filter) ' km'])
        xlabel('Time')
        ylabel('Hs (m)')
        xlim([datetime(1979,01,01) datetime(2010,12,31)])
        fontsize(fig_timehistory,12,"points")

        %% plot the seasonality
        fig_seasonality = figure;
        fig_seasonality.Position = [100 100 800 300];
        s3 = scatter(month(TC_peaks.Time),TC_peaks.Hs,Marker="x",DisplayName="TC peaks",SizeData=60,MarkerEdgeColor="blue",LineWidth=2);
        hold on
        s4 = scatter(month(extremes_TT.Time),extremes_TT.Hs,Marker="o",DisplayName="nonTC peaks",SizeData=60,MarkerEdgeColor="red",MarkerFaceColor='red');
        xticks(1:12); % Set x-ticks to represent each month
        xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'}); % Label x-ticks
        xlabel('Month');
        ylabel('Hs (m)');
        % title([sites.Location{loc} ' TC and nonTC peaks month R = ' num2str(distance_filter) ' km'])
        xlim([0 13])
        grid on
        fontsize(fig_seasonality,16,"points")
        legend(Location="best",FontSize=12)


        %% export all figure
        exportgraphics(fig_timehistory,sprintf('%s Hs timehistory %0.0f radius all hurdat.png',site_desc,distance_filter),Resolution=450)
        exportgraphics(fig_seasonality,sprintf('%s seasonality %0.0f radius all hurdat.png',site_desc, distance_filter),Resolution=450)

        %% save the neeeded parameters
        IID_NREL_all_HURDAT(loc).name         = site_desc;
        IID_NREL_all_HURDAT(loc).Hs_all       = Hs_data;
        IID_NREL_all_HURDAT(loc).TC_data      = TC_data_filtered;
        IID_NREL_all_HURDAT(loc).MIS_TC       = TC_peaks;
        IID_NREL_all_HURDAT(loc).nonTC_data   = nonTC_data;
        IID_NREL_all_HURDAT(loc).nonTC_AM     = extremes_TT;
        close all
    end
    alldist_all_hurdat(dista).result = IID_NREL_all_HURDAT;
    alldist_all_hurdat(dista).searchdist = dist_array(dista);
end

save("alldist_all_hurdat","alldist_all_hurdat")