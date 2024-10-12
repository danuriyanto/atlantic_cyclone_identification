clc; clear; close all
addpath cyclone_id_functions/
load("alldist_hu_only.mat")
load("stats_method_fitting_HU_only.mat")
IID_NREL = alldist_hu_only(5).result;

nonTC_params = array2table([0 0 0],VariableNames={'k' 'sigma' 'mu'});
TC_params    = nonTC_params;

for loc = 1:12
    Hs_all = IID_NREL(loc).Hs_all;
    th = prctile(Hs_all.Hs,90);

    %% get the nonTC duration
    AM = IID_NREL.nonTC_AM;
    [duration_nonTC,~] = calc_duration(AM,Hs_all);
    duration_nonTC = rmmissing(duration_nonTC);

    %% get TC duration
    MIS = IID_NREL(loc).MIS_TC;
    [duration_TC,~] = calc_duration(MIS,Hs_all);
    duration_TC = rmmissing(duration_TC);

    %% fit into distribution
    [dist_nonTC] = dist_duration(duration_nonTC);
    R2_nonTC = calc_R2((duration_nonTC),dist_nonTC)
    [dist_TC] = dist_duration(duration_TC);
    R2_TC = calc_R2((duration_TC),dist_TC)

    %% parameters
    nonTC_params{loc,:} = dist_nonTC.ParameterValues;
    TC_params{loc,:} = dist_TC.ParameterValues;

    %% PoE
    T = [1:1:200];
    PoE_TC    = 1 - cdf(dist_TC, T);
    PoE_nonTC = 1 - cdf(dist_nonTC,T); 

    f = figure;
    semilogx(PoE_TC,T,LineWidth=1.2,LineStyle="--",Color=[0 0 0], ...
        DisplayName='TC')
    hold on
    semilogx(PoE_nonTC,T,LineWidth=1.2,LineStyle="-",Color=[0.4 0.4 0.4], ...
        DisplayName='nonTC')
    yticks(0:24:24*8)
    legend(Location="northeast")
    grid on
    xlabel('PoE')
    ylabel('duration (hours)')
    title(sprintf('%s',alldist_hu_only(1).result(loc).name))
    fontsize(f,14,"points")
    exportgraphics(f,['PoE_' alldist_hu_only(1).result(loc).name '.png'],Resolution=450)
end
%% functions
function [dist_result] = dist_duration(duration_vector)
dist_result = fitdist(duration_vector,"GeneralizedExtremeValue");
xx = min(duration_vector)-20:5:max(duration_vector)+20;
pdf_values = pdf(dist_result,1:200);

% figure
% histogram(duration_vector,10,Normalization="pdf")
% hold on
% plot(1:200,pdf_values)
% 
% [aa,bb] = ecdf(duration_vector);
% cdf_val = cdf(dist_result,xx);
% 
% figure
% stairs(bb,aa)
% hold on
% plot(xx,cdf_val)
end