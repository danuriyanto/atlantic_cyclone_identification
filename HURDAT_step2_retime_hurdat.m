clc; clear; close all;
load("hurdat2_track_1979_2010.mat")
addpath cyclone_id_functions/

hur_id = unique(hurdat2_track_1979_2010.hurricane_yr_name);
hurdat2_retimed = cell([numel(hur_id) 1]);
for hur=1:numel(hur_id)
    % load the hurricane TimeTable unique to their name
    TT = hurdat2_track_1979_2010(hurdat2_track_1979_2010.hurricane_yr_name == hur_id(hur),:);

    %% cleanup the hurricane
    TT = sortrows(TT,"Time");   %sort the time
    TT = unique(TT);            %cleanup the duplicate values
    
    % Find Rows with Duplicate Times and Different Data
    dupTimes = sort(TT.Time);
    tf = (diff(dupTimes) == 0);
    dupTimes = dupTimes(tf);
    dupTimes = unique(dupTimes);

    if ~isempty(dupTimes)
        continue
    end
    % 
    % % find the first value of the duplicates
    % uniqueTimes = unique(TT.Time);
    % TT = retime(TT,uniqueTimes,"mean");
    % 
    Time = (TT.Time(1):hours(3):TT.Time(end))';

    % Interpolate numeric data
    lon_interp = interp1(TT.Time, TT.lon, Time, 'spline');
    lat_interp = interp1(TT.Time, TT.lat, Time, 'spline');
    u10_interp = interp1(TT.Time, TT.u10_max_mps, Time, 'spline');

    % Duplicate string data for new time vector
    storm_ID_new = strings(size(Time));
    status_new = strings(size(Time));
    for i = 1:length(Time)
        % Find the original time point closest to the current new time point
        [~, idx] = min(abs(TT.Time - Time(i)));
        storm_ID_new(i) = TT.hurricane_yr_name(idx);
        status_new(i) = TT.status{idx};
    end

    % Create new timetable
    TT_new = timetable(Time, storm_ID_new, status_new, lon_interp, lat_interp, u10_interp);
    TT_new.Properties.VariableNames = {'storm_ID', 'status', 'lon', 'lat', 'u10'};

    hurdat2_retimed{hur,:} = TT_new;
end

hurdat2_retimed = vertcat(hurdat2_retimed{:});
save("hurdat2_retimed_1979_2010_2","hurdat2_retimed")
% hold on
% text(TT_new.lat,TT_new.lon,TT_new.status)