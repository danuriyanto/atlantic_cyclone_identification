function extreme_vector = get_r_yearlymax(timehistory,r)
extreme_vector = [];
minyear = min(year(timehistory.Time));
maxyear = max(year(timehistory.Time));
dt = hours(72);
for yr=minyear:maxyear
    startdate = datetime(sprintf('1/1/%d 00:00:00',yr),"InputFormat","MM/dd/uuuu HH:mm:ss");
    enddate = datetime(sprintf('12/31/%d 23:59:00',yr),"InputFormat","MM/dd/uuuu HH:mm:ss");
    tr = timerange(startdate,enddate);
    extremes_th_at_year_yr = timehistory(tr,:); %yearly value of non-retimed Hs
    selection = extremes_th_at_year_yr;
    for i=1:r
        [r_extreme,idhs] = max(selection.Hs); %take the maximum of retimed Hs
        datemax = selection(idhs,:); %find the date of retimed Hs
        daterange = datemax.Time-dt/2:hours(1):datemax.Time+dt/2;
        selection(daterange,"Hs") = table(zeros(length(daterange),1));
        selected_extreme(i,:) = r_extreme;
        id(i,:) = idhs;
    end
    extreme_vector = [extreme_vector; selected_extreme];
    extreme_vector = sort(extreme_vector,"ascend");
end

end
