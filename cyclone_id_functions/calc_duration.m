function [duration_vector,mean_duration] = calc_duration(storm_TT,Hs_all_TT)
% CALC_DURATION Calculate the duration of extreme events on
% written by: Danu Riyanto, 2024-08-23
% Inputs:
%   TT      - Timetable containing specific times to evaluate in the Time field
%   Hs_all  - Timetable containing all times and corresponding significant wave heights (Hs)
%
% Outputs:
%   duration_vector - Vector of durations (in hours) for each event above the 90th percentile threshold
%   mean_duration   - Mean duration of all events in the duration_vector

% Initialize the duration vector
duration_vector = NaN([size(storm_TT,1) 1]);

% Loop through each year (or entry) in the TT table
for yr = 1:size(storm_TT,1)
    % Calculate the 90th percentile threshold for significant wave height (Hs)
    th = prctile(Hs_all_TT.Hs,90);
    
    % Find the index of the maximum Hs value at the specified time in TT
    id_max = find(Hs_all_TT.Time == storm_TT.Time(yr));
    
    % Skip the iteration if the Hs at the specified time is below the threshold
    if Hs_all_TT.Hs(id_max) < th
        continue
    end
    
    % Initialize the forward and backward search indices and peak Hs values
    id_forward = id_max;
    id_backward = id_max;
    peakforward = Hs_all_TT.Hs(id_max);
    peakbackward = peakforward;
    
    % Initialize arrays to store the indices of the event duration
    id_forward_all = [];
    id_backward_all = [];
    
    % Search for the duration of the event by expanding forward and backward
    while peakforward > th && peakbackward > th
        % Move backward until Hs drops below the threshold
        id_backward = id_backward - 1;
        peakbackward = Hs_all_TT.Hs(id_backward);
        id_backward_all = [id_backward_all; id_backward];

        % Move forward until Hs drops below the threshold
        id_forward = id_forward + 1;
        peakforward = Hs_all_TT.Hs(id_forward);
        id_forward_all = [id_forward_all; id_forward];
    end

    % Combine forward and backward indices, sort them, and calculate the duration
    id_storm = sort([id_backward_all; id_forward_all], "ascend");
    nonTC_TT = Hs_all_TT.Time(id_storm);
    duration_vector(yr,:) = hours(nonTC_TT(end) - nonTC_TT(1));
end

% Calculate the mean duration of all identified events
mean_duration = mean(duration_vector);