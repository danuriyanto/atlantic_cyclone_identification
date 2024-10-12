function [MRP,MRP_selection,distparameter] = annual_maxima_rlargest(extremes_timehistory, r, grouptitle,plotflag,dist_type,RP_selection)

%predefined input within the function
maxyear = year(max(extremes_timehistory.Time));
minyear = year(min(extremes_timehistory.Time));
dt = hours(24*3);
extreme_vectors = [];
plotselections = timetable();

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
        plotselection = extremes_th_at_year_yr(id,:);
    end
    plotselections = [plotselections; plotselection];
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

[MRP,MRP_selection,distparameter] = calc_annual_maxima_rlargest(...
    extreme_vectors,r,grouptitle,plotflag,dist_type,RP_selection,maxyear,minyear);

end

function [MRP,MRP_selection,distparameter] = calc_annual_maxima_rlargest(...
    extremes_vector,r,grouptitle,plotflag,dist_type,RP_selection,maxyear,minyear)

%fitting into distributions
if dist_type == "GEV"
    distparameter     = fitdist(extremes_vector,"GeneralizedExtremeValue");
elseif dist_type == "Weibull"
    distparameter     = fitdist(extremes_vector,"Weibull");
elseif dist_type == "Gumbel"
    distparameter     = fitdist(extremes_vector,"ExtremeValue");
elseif dist_type == "Exponential"
    distparameter     = fitdist(extremes_vector,"Exponential");
elseif dist_type == "GEV_limit_k"
    [~,~,~,distparameter] = makedist_gev_limit(extremes_vector);
else
    error('no such distribution goblog!!')
end


%calculating the return value
maxT    = 100000;
T       = linspace(1,maxT,maxT)';
RP      = 1-(1./(r*T));
MRI     = icdf(distparameter,RP);
MRP     = [T MRI];
MRP_selection = MRP(RP_selection,:);

%calculating quantiles
n       =length(extremes_vector);
qq      = sort(extremes_vector);
p       = ((1:n) - 0.5) ./ n;
q       = icdf(distparameter,p)';
R2      = 1 - sum((qq - q).^2)/sum((qq - mean(qq)).^2);

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
    fig.Position = [100 100 900 700];
    fig.Color = 'white';
    subplot(2,2,1)
    histogram(extremes_vector,20,Normalization="pdf",FaceColor='b',FaceAlpha=.2,EdgeColor='none')
    hold on
    plot(x,pdf(distparameter,x),"b",'LineWidth',.9)
    grid on
    xlabel(grouptitle)
    ylabel('pdf')
    legend('histogram',dist_type,Location='northeast')
    title('histogram vs pdf')

    %plotting empirical cdf vs theoretical cdf
    [f1,x1]=ecdf(extremes_vector);
    subplot(2,2,2)
    stairs(x1,f1,'k','LineWidth',.9)
    hold on
    plot(x,cdf(distparameter,x),"b",'LineWidth',.9)
    grid on
    xlabel(grouptitle)
    ylabel('cdf')
    legend('empirical cdf',dist_type,Location='southeast')
    title('empirical vs theoretical cdf')

    %plotting return periods
    [C,~,ic] = unique(extremes_vector,'sorted'); % ic are ranks from lowest to highest ; C are unique values
    rr=(1+max(ic)-ic);  % r: rank (highest receives 1; lowest receives length(C); tied values receive same rank)
    nn = maxyear+1-minyear;
    P = (rr)./(nn+1);
    RP = 1./P;
    RP = sort(RP,'ascend');

    subplot(2,2,3)
    semilogx(T,MRI,"b",'LineWidth',.9)
    xlim([0 100000])
    box on
    hold on
    grid on
    semilogx(RP,sort(extremes_vector,"ascend"),'-ok',MarkerSize=3,MarkerFaceColor='black',LineWidth=.4)
    legend(dist_type,'actual return period',Location='northwest')
    title('return period')
    xlabel('return period')
    ylabel('Hs (m)')
    ylim([0 ceil(max(MRI))+10])

    %plotting Q-Q
    subplot(2,2,4)
    box on
    hold on
    grid on
    plot(qq,q,'ob',MarkerSize=3,MarkerFaceColor='b')
    plot(qq,qq,'k',LineWidth=1)
    title('Q-Q plot')
    xlabel('theoretical')
    ylabel('observed')
    xlim([min(qq)-1 max(qq)+1])
    ylim([min(qq)-1 max(qq)+1])
    xticks(floor(min(qq)):1:(ceil(max(qq))+2))
    yticks(floor(min(qq)):1:(ceil(max(qq))+2))
    R2label = sprintf('%s (R^2=%.3f)',dist_type, R2);
    legend(R2label,Location='southeast')
    fontsize(fig,14,"points")
    if plotflag == 'y' && ~isempty(grouptitle)
        if     r==1
            exportgraphics(fig,join([grouptitle '-' dist_type '-AM.png']),Resolution=300)
        elseif r>=2
            exportgraphics(fig,join([grouptitle '-' dist_type '-R2.png']),Resolution=300)
        end
    end
end


end