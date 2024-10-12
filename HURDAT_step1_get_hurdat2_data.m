clear
clc
addpath cyclone_id_functions/
file = ('data/hurdat2-1851-2023-051124.txt');
hurricane_id = readtable("data/hurricane_ID_1979-2010.csv",ReadVariableNames=false);
yr = zeros(size(hurricane_id));
table_header = {'hurricane_yr_name','status','lon','lat','u10_max_mps'};
for i=1:size(hurricane_id,1)
    i
    [t{i}, ~, status{i}, lon{i}, lat{i}, kt_max{i}, ~, ~, ~, ~, ~, ~, ~,...
        ~, ~, ~, ~, ~, ~, ~, hurricane_name{i}] ...
        = readHurdat2(file, hurricane_id{i,1});
    kt_max{i} = kt_max{i}*0.51444;
    yr(i) = year(t{i}(1));
    hur_yr_name{i} = [num2str(yr(i)) '_' hurricane_name{i}];
    hur_yr_name_array{i} = repmat(convertCharsToStrings(hur_yr_name{i}),size(kt_max{i}));
    hur = timetable(t{i},hur_yr_name_array{i}, status{i},lon{i}, lat{i}, kt_max{i}, VariableNames=table_header);
    % hur.hurricane_yr_name = categorical(hur.hurricane_yr_name);
    hurdat2_track_1979_2010{i,:} = hur; 
end

hurdat2_track_1979_2010 = vertcat(hurdat2_track_1979_2010{:});
save("hurdat2_track_1979_2010.mat","hurdat2_track_1979_2010");