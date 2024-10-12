function stats_method_result = stats_method(data, observation_year, distribution_obj, RP, sitename, storm_type,plotdiagnostics)
% statistical method to check the MRP

%% calculations
% empirical CDF
n=length(data);
xx = sort(data);
ff = ((1:n)-0.5)' ./ n;
x_data = floor(min(data)):0.01:ceil(max(data));

%calculating quantiles
p = ((1:n) - 0.5) ./ n;
q_LS = icdf(distribution_obj,p)';

% Empirical percentiles (CDF values)
epercentiles      = (1:length(data)) / length(data);

% Theoretical percentiles from the fitted GPD object
tpercentiles_LS   = cdf(distribution_obj, data);

% calculating MRP
maxT                 = 100000;                          % calculate max desired return value
T                    = linspace(0.01,maxT,maxT)';       % calculate array of years
aT                   = [0.01:0.01:3]';
T                    = sort(unique([T; aT]),"ascend");
lam                  = size(data,1)/(observation_year); % arrival rate
PoE                  = 1-(1./(lam*T));                  % probability of exceedance
MRI                  = icdf(distribution_obj,PoE);      % MRI from fitted distribution
MRP_interp           = round(interp1(T,MRI,RP,"linear"),2);
MRI_requested        = [RP MRP_interp];
MRP                  = [T MRI];

% calculating empirical MRP
[RP,MRI_empirical] =  plotMRI(data,observation_year);

empirical = [RP MRI_empirical];

%% plottings
% plotting CDF
if plotdiagnostics == 'y'
    fig_diagnostics = figure(Visible="off");
    fig_diagnostics.Position = [200 300 1120 840];
    subplot(221)
    stairs(xx,ff, LineWidth=1, Color='black',DisplayName='ecdf TC')
    hold on
    plot(x_data,cdf(distribution_obj,x_data),LineWidth=1, Color='b',DisplayName='fitted CDF')
    legend(Location="southeast")
    xlabel('Hs (m)')
    ylabel('CDF')

    % plotting return period
    subplot(222)
    semilogx(T,MRI,'--b', LineWidth=1,DisplayName='MRP')
    hold on
    semilogx(RP,MRI_empirical,LineStyle="--",Marker="o",Color='k',DisplayName='empirical MRI')
    legend(Location="southeast")
    xlabel('return period')
    ylabel('return level - Hs (m)')
    xlim([1 100000])
    xticks([0 1e1 1e2 1e3 1e4 1e5])

    % plotting quantiles and R2 quantiles
    R2 = calc_R2(data, distribution_obj);
    subplot(223)
    plot(data,data,'--k',DisplayName='reference', LineWidth=1)
    hold on
    plot(q_LS,data,'ob',DisplayName=sprintf('QQ, R^2 = %0.03f',R2), LineWidth=1)
    xlabel('theoretical')
    ylabel('observed')
    legend(Location="southeast")

    % plotting percentiles
    subplot(224)
    plot(tpercentiles_LS,tpercentiles_LS,'--k',DisplayName='reference',LineWidth=1)
    hold on
    plot(tpercentiles_LS,epercentiles,'ob',DisplayName='PP', LineWidth=1)
    xlabel('theoretical')
    ylabel('observed')
    legend(Location="southeast")
    sgtitle([sitename ' ' storm_type ' diagnostic'])
    fontsize(fig_diagnostics,16,"points")
    exportgraphics(fig_diagnostics,[sitename '_' storm_type '.png'])
end

% save into struct 
stats_method_result = struct();
stats_method_result.MRP_fitted     = MRP;
stats_method_result.MRP_empirical  = empirical;
stats_method_result.MRI_requested  = MRI_requested;
stats_method_result.empirical_cdf  = [xx ff];
% stats_method_result.R2_quantiles   = R2;
stats_method_result.dist_type      = storm_type;
stats_method_result.dist_object    = distribution_obj;
stats_method_result.ks_test        = kstest(data,"CDF",[data cdf(distribution_obj, data)]);
stats_method_result.ad_test        = adtest(data,Distribution=distribution_obj);

end