function [PoT_dataset,thres] = get_PoT_data(timehistory,threshold_percentile,dt)

% retiming the data based on assumed storm duration
timehistory_retimed = retime(timehistory,"regular","max",TimeStep=dt);

% generating the Peak over Thershold Dataset
[peaksval, ~] = findpeaks(timehistory_retimed{:,1});
thres = prctile(peaksval,threshold_percentile);
[~, peaksloc] = findpeaks(timehistory.Hs,timehistory.Time,MinPeakHeight=thres,MinPeakDistance=dt);
th_line = ones(size(timehistory,1),1)*thres;
th_plot = timetable(timehistory.Time,th_line);
PoT_dataset = timehistory(peaksloc,1);
