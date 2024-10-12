clc; clear; close all;
addpath cyclone_id_functions/
load('data/hurdat2_retimed_1979_2010.mat')
sites = readtable("data/12_sites_final.csv");

% load Hs and local u10 from 12 sites
Hs_data = load("data/EVA_METHOD_COMPARISON_2024/Hs_NREL.mat");
Hs_data = Hs_data.Hs_NREL;
Hs_data = retime(Hs_data,"hourly","linear");
load("data/V10_NCEP.mat");

hur_name = unique(hurdat2_retimed.storm_ID);
wgs = wgs84Ellipsoid("kilometer");

for loc = 1:size(sites,1)
    %% get all the passing TC
    clearvars -except Hs_data Vdata hur_name wgs dt sites loc hurdat2_retimed

    TC_data = cell([length(hur_name), 1]);
    lon     = sites.Longitude(loc);
    lat     = sites.Latitude(loc);
    Hs      = Hs_data(:,loc);
    null_hs = Hs{:,1} <=0;
    Hs      = Hs(~null_hs,:);

    for n=1:length(hur_name)
        hur_at_n        = hurdat2_retimed...
            (hurdat2_retimed.storm_ID==hur_name(n),:);
        dist            = distance(lat,lon,hur_at_n.lat,hur_at_n.lon,wgs);
        if any(dist<=1000)
            sprintf('found at n = %0.0f',n)
            % retime the hurricane timeseries and slice the NREL Hs data
            hur = hur_at_n(:,2:end);
            hur = sortrows(hur,"Time");
            uniquetimes = unique(hur.Time);
            null_val    = hur.u10<=0;
            hur         = hur(~null_val,:);
            Hs_hur      = Hs{hur.Time,:};
            V_local_hur = Vdata{hur.Time,loc};

            % add the Hs and categorical name to retimed hurricane dataset
            dstnce = distance(lat,lon,hur.lat,hur.lon,wgs);
            hur = addvars(hur,...
                repmat(string(hur_at_n.storm_ID(1)),[size(hur,1),1]),...
                Before='lon',NewVariableNames='storm_ID');
            hur = addvars(hur, ...
                dstnce,NewVariableNames='distance_from_position');
            hur = addvars(hur,Hs_hur,NewVariableNames='Hs');
            hur = addvars(hur,V_local_hur,NewVariableNames='u10_local');
            TC_data{n} = hur;
        end
    end
    TC_data = vertcat(TC_data{:});

    %% get the duplicate tracks
    maxHs = groupsummary(TC_data,"storm_ID","max","Hs");
    maxHs = sortrows(maxHs,"max_Hs","descend");
    maxHs.max_Hs = round(maxHs.max_Hs,3);

    for n = 1:height(maxHs)
        % get the storm date when maxHs occured.
        storm_at_ID_n = TC_data(TC_data.storm_ID == maxHs.storm_ID(n),:);
        [~,idmax] = max(storm_at_ID_n.Hs);
        storm_date(n,:) = storm_at_ID_n.Time(idmax);

        % get the minimum distance when maxHs occured to check which storm
        % caused the maximum Hs, because it is possible that we have
        % multiple storms in one instances.
        storm_distance(n,:) = storm_at_ID_n.distance_from_position(idmax);
    end

    % add the data and distance of the closest storm eye to the turbine pos
    maxHs = addvars(maxHs,storm_date);
    maxHs = addvars(maxHs,storm_distance);

    % get the indices of the duplicate storms
    [C, ia, ic] = unique(maxHs.storm_date); % get the index of unique storm date
    dupl_Hs_id = setdiff(1:size(maxHs,1),ia(accumarray(ic,1)<=1))'; % get the duplicate index
    maxHs_multiple_storms = maxHs(dupl_Hs_id,:); % group the storm metadata based on their duplicates

    % get the minimum distance and disregard the rest
    [b,bc,bg] = groupsummary(maxHs_multiple_storms.storm_distance,maxHs_multiple_storms.storm_date,"min","IncludeEmptyGroups",true);
    for dup = 1:length(bc)
        dup_date = bc(dup);
        dup_occurence = maxHs_multiple_storms(maxHs_multiple_storms.storm_date==dup_date,:);
        [qq,rr] = min(dup_occurence.storm_distance);
        dup_occurence(rr,:) = [];
        st_ID_to_remove = dup_occurence.storm_ID;
        for nn = 1:length(st_ID_to_remove)
            NaNs = NaN(size(TC_data(TC_data.storm_ID==st_ID_to_remove(nn),"Hs")));
            TC_data.Hs(TC_data.storm_ID==st_ID_to_remove(nn)) = NaNs;
        end
    end

    save(join(['TC_' sites.Location{loc}]),"TC_data")
end

