function [extremes_TT,extreme_vectors,distparameter] = rlargest_CDF(extremes_timehistory, r, grouptitle,plotflag,dist_type,RP_selection)

%predefined input within the function
maxyear = year(max(extremes_timehistory.Time));
minyear = year(min(extremes_timehistory.Time));
dt = hours(24*3);
extreme_vectors = [];
extremes_TT = timetable();

%looping for each year to get the r yearly maxima
for yr=minyear:maxyear
    startdate = datetime(sprintf('1/1/%d 00:00:00',yr),"InputFormat","MM/dd/uuuu HH:mm:ss");
    enddate = datetime(sprintf('12/31/%d 23:59:00',yr),"InputFormat","MM/dd/uuuu HH:mm:ss");
    tr = timerange(startdate,enddate);
    extremes_th_at_year_yr = extremes_timehistory(tr,:); %yearly value of non-retimed Hs
    selection = extremes_th_at_year_yr;
    for i=1:r
        [r_extreme,idhs] = max(selection.Hs); %take the maximum of retimed Hs
        datemax = selection(idhs,:); %find the date of retimed Hs
        daterange = datemax.Time-dt/2:hours(1):datemax.Time+dt/2;
        selection(daterange,"Hs") = table(zeros(length(daterange),1));
        selected_extreme(i,:) = r_extreme;
        id(i,:) = idhs;
        plotselection = extremes_th_at_year_yr(idhs,:);
    end
    extremes_TT = [extremes_TT; plotselection];
    extreme_vectors = [extreme_vectors; selected_extreme];
    %plot all the occurences in all location (commented out bcs already did)
    % figr = figure(Visible="on");
    % figr.Position = [400 400 800 450];
    % hold on
    % plot(extremes_th_at_year_yr.Time,extremes_th_at_year_yr{:,1},'k')
    % plot(plotselection.Time,plotselection{:,:},'xr')
    % grid on
    % xlabel('date')
    % ylabel('Hs (m)')
end

% fig1 = figure;
% fig1.Position = [200 200 1200 400];
% plot(year(plotselections.Time,"iso"),plotselections{:,1},MarkerSize=20,LineStyle="none",Marker=".")
% xticks(1979:1:2010)
% xlim([1979 2010])
% ylim([2 14])
% grid on
% box on
% xlabel('year')
% ylabel('Hs (m)')
% fontsize(fig1,16,"points")

[distparameter] = calc_annual_maxima_rlargest(extreme_vectors,r,grouptitle,plotflag,dist_type,RP_selection);

end

function [distparameter] = calc_annual_maxima_rlargest(extremes_vector,r,grouptitle,plotflag,dist_type,RP_selection)

%fitting into distributions
if dist_type == "GEV"
    distparameter     = fitdist(extremes_vector,"GeneralizedExtremeValue");
elseif dist_type == "Weibull"
    distparameter     = fitdist(extremes_vector,"Weibull");
elseif dist_type == "Gumbel"
    distparameter     = fitdist(extremes_vector,"ExtremeValue");
elseif dist_type == "Exponential"
    distparameter     = fitdist(extremes_vector,"Exponential");
else
    error('no such distribution!!')
end

if plotflag == "y"
    %plot histogram vs pdf plot
    x   = linspace(0,20,100);
    fig = figure;
    if isempty(grouptitle)
        sgtitle(' ');
    elseif r == 1
        sgtitle(join([grouptitle 'annual maxima fitted into' dist_type]))
    elseif r >= 2
        sgtitle(join([grouptitle 'r =' num2str(r) 'fitted into' dist_type]))
    end
    fig.Position = [100 100 800 400];
    fig.Color = 'white';

    %plotting empirical cdf vs theoretical cdf
    [f1,x1]=ecdf(extremes_vector);
    subplot(2,3,[1 2 4 5])
    stairs(x1,f1,'k','LineWidth',.9)
    hold on
    plot(x,cdf(distparameter,x),"b",'LineWidth',.9)
    grid on
    xlabel('Hs (m)')
    ylabel('cdf')
    % legend('empirical cdf',dist_type,Location='southeast')
    title('empirical vs theoretical cdf')

    subplot(2,3, [3 6])
    stairs(x1,f1,'k','LineWidth',.9)
    hold on
    plot(x,cdf(distparameter,x),"b",'LineWidth',.9)
    grid on
    xlabel('Hs (m)')
    legend('empirical cdf',dist_type,Location='southoutside')
    title('upper tail')
    xlim([prctile(x1,80) prctile(x1,100)*1.5])
    ylim([prctile(f1,90) prctile(f1,100)])

    fontsize(fig,16,"points")
    if plotflag == 'y' && ~isempty(grouptitle)
        if     r==1
            exportgraphics(fig,join([grouptitle '-' dist_type '-AM.png']),Resolution=300)
        elseif r>=2
            exportgraphics(fig,join([grouptitle '-' dist_type '-R2.png']),Resolution=300)
        end
    end

end
end