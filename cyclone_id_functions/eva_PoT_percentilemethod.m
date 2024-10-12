function [MRP,MRP_selection,distribution_pararameter,R2] = ...
    eva_PoT_percentilemethod(extremes, dist_type, pct, RP_selection, plot_show_flag,plot_save_flag)
%To calculate the extreme value of a timeseries using Peak over Threshold
%(PoT) using percentile method. Written by Raditya Danu Riyanto, 6/18/2023
% The peaks are chosen after declustering, with storm duration = 3 days. 
% The inputs are:
% 1. timeseries_table: Timeseries of the hourly observed value (Matlab
%                      Timetable format) 
% 2. pct             : threshold percentile value (int/float)
% 3. T_selection     : Return period of interest(array 1-1000)
% 4. plot_show_save  : flag to show the diagnostic plot ("y"/"n")

%%%%%%%%%%%%%%%%%% predefined values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titl = extremes.Properties.VariableNames;
dt = hours(3);

%%%%%%%%%%%%%%%%%% Generating the data above threshold %%%%%%%%%%%%%%%%%%%%

% retiming the data based on assumed storm duration
extremes_retimed = retime(extremes,"regular","max",TimeStep=dt);

% generating the Peak over Thershold Dataset
[~, peaksloc] = findpeaks(extremes_retimed{:,1},"MinPeakHeight",pct);
th_line = ones(size(extremes_retimed,1),1)*pct;
th_plot = timetable(extremes_retimed.Time,th_line);
PoT_dataset = extremes_retimed(peaksloc,1);


%plotting the Peak over Thershold Dataset
if plot_show_flag=="y"
    fig1 = figure(Visible="on");
else
    fig1 = figure(Visible="off");
end
fig1.Position = [300 200 1100 400];
hold on
plot(extremes_retimed.Time,extremes_retimed{:,1},LineWidth=0.3)
plot(th_plot.Time,th_line,'--',LineWidth=1.3,Color='magenta')
plot(PoT_dataset.Time,PoT_dataset{:,1},'.r')
ylim([0 max(extremes{:,1})+2])
grid on
xlabel('year')
ylabel(titl)
title(join([string(titl) sprintf("PoT based on 90th percentile threshold = %.1f Dataset",pct)]))
if plot_save_flag == "y"
    exportgraphics(fig1,join([string(titl) sprintf("PoT based on 90th th percentile threshold = %.1f Dataset.png", pct)]),Resolution=600)
end

%%%%%%%%%% Fitting into distributions and calculate return periods %%%%%%%%
if dist_type == "GPD" || dist_type == "gpd" || dist_type == "generalized pareto"
    distribution_pararameter = fitdist(PoT_dataset{:,1},"GeneralizedPareto",Theta=pct_val);
elseif dist_type == "weibull" || dist_type == "Weibull" || dist_type == "wbl" || dist_type == "WBL"
    distribution_pararameter = fitdist(PoT_dataset{:,1},"Weibull");
elseif dist_type == "Exponential" || dist_type == "exponential" || dist_type == "exp" || dist_type == "EXP"
    distribution_pararameter = fitdist(PoT_dataset{:,1},"Exponential");
elseif dist_type == "GPD_limit"
    [~,~,~,distribution_pararameter] = makedist_gpd_limit(PoT_dataset{:,1});
end


%calculating the return value
maxT  = 100000;
T     = linspace(1,maxT,maxT)';
obsv_duration =  max(year(PoT_dataset.Time))+1 - min(year(PoT_dataset.Time));
lam   = size(PoT_dataset,1)/(obsv_duration);
RP    = 1-(1./(lam*T));
MRI   = icdf(distribution_pararameter,RP);

%calculating quantiles
n=length(PoT_dataset{:,1});
qq = sort(PoT_dataset{:,1});
p = ((1:n) - 0.5) ./ n;
q_dist = icdf(distribution_pararameter,p)';
R2 = 1 - sum((qq - q_dist).^2)/sum((qq - mean(qq)).^2);


%%%%%%%%%%%%%%%%%%%%%%%% constructing diagnostics plot %%%%%%%%%%%%%%%%%%%%
%plot histogram vs pdf plot
x = linspace(0,max(extremes{:,1})+4,100);
if plot_show_flag=="y"
    fig2 = figure(Visible="on");
else
    fig2 = figure(Visible="off");
end

sgtitle(join([string(titl) sprintf("PoT based on 90th percentile threshold = %.1f m fitted into %s Diagnostic",pct, dist_type)]))
fig2.Position = [100 100 900 700];
fig2.Color = 'white';

subplot(2,2,1)
histogram(PoT_dataset{:,1},20,Normalization="pdf",FaceColor='b',FaceAlpha=.2,EdgeColor='none')
hold on
plot(x,pdf(distribution_pararameter,x),"r",'LineWidth',.9)
grid on
xlabel(titl)
ylabel('pdf')
legend('histogram',dist_type,Location='northeast')
title('histogram vs pdf')

%plotting empirical cdf vs theoretical cdf
[f1,x1]=ecdf(PoT_dataset{:,1});
subplot(2,2,2)
stairs(x1,f1,'k','LineWidth',.9)
hold on
plot(x,cdf(distribution_pararameter,x),"r",'LineWidth',.9)
grid on
xlabel(titl)
ylabel('cdf')
legend('empirical cdf',dist_type,Location='southeast')
title('empirical vs theoretical cdf')

%plotting return periods
[C,~,ic] = unique(PoT_dataset{:,1},'sorted'); % ic are ranks from lowest to highest ; C are unique values
rr=(1+max(ic)-ic);  % r: rank (highest receives 1; lowest receives length(C); tied values receive same rank)
nn = obsv_duration;
P = (rr)./(nn+1);
RP = 1./P;
RP = sort(RP,'ascend');

subplot(2,2,3)
semilogx(T,MRI,"r",'LineWidth',.9)
box on
hold on
grid on
plot(RP,sort(PoT_dataset{:,1},"ascend"),'-ok',MarkerSize=3,MarkerFaceColor='black',LineWidth=.4)
legend(dist_type,'actual return period',Location='best')
title('return period')
xlabel('return period')
ylabel(titl)

%plotting Q-Q
subplot(2,2,4)
box on
hold on
grid on
plot(qq,q_dist,'or',MarkerSize=3,MarkerFaceColor='r')
plot(qq,qq,'k',LineWidth=1)
title('Q-Q plot')
xlabel('theoretical')
ylabel('observed')
xlim([min(qq)-1 max(qq)+1])
ylim([min(qq)-1 max(qq)+1])
xticks(floor(min(qq)):1:(ceil(max(qq))+2))
yticks(floor(min(qq)):1:(ceil(max(qq))+2))
l_dist = sprintf('GPD(R^2=%.3f)',R2);
legend(l_dist,Location='best')

if plot_save_flag == "y"
    exportgraphics(fig2,join([string(titl) sprintf("PoT based on 90 th percentile threshold = %.1f m fitted into %s Diagnostic.png",pct,dist_type)]),Resolution=600)
end

%%%%%%%%%%%%%%%%%%% creating return value table and plot %%%%%%%%%%%%%%%%%%
MRP = [T MRI];
MRP_selection = MRP(RP_selection,:);

end