function empiricalResult = empiricalMRP(lambda_TC,lambda_nonTC,dist_TC,dist_nonTC,RP)
% #1: directly draw 100,000 data and plot the MRI
rng(12345)
yr = 20000000;

%% draw the random number along the target MRP year for nonTC
arr_nonTC = (lambda_nonTC*yr);
point_drawn_nonTC = random(dist_nonTC, [arr_nonTC 1]);



%% draw the random number along the target MRP year for TC
arr_TC = int32(lambda_TC*yr);
point_drawn_TC = random(dist_TC,[arr_TC 1]);

% Initialize a 100,000xn array with zeros
R = zeros(yr, randi([2 20],1));
num_R = numel(R); % get the number of element
randp = randperm(num_R,arr_TC); % random permutation number 
tall_array = zeros([num_R,1]);  % generate a tall array with size num_R
tall_array(randp) = point_drawn_TC; % reshape to the original R
points_TC_all = reshape(tall_array,size(R));

%% combine the TC
points_combined = max([point_drawn_nonTC points_TC_all],[],2);
point_all = sort(unique(points_combined), "ascend");
[T,MRI] = plotMRI(point_all,yr);
MRI_request = round(interp1(T,MRI,RP,"linear"),2);
MRI_table = table(T,MRI);
empiricalResult = struct;
empiricalResult.RP_request = table(RP,MRI_request);
empiricalResult.MRI_table  = MRI_table;

%% functions local
    function [RP,MRI] =  plotMRI(maximas,observation_year)
        [C,~,ic] = unique(maximas,'sorted'); % ic are ranks from lowest to highest ; C are unique values
        rr=(1+max(ic)-ic);  % r: rank (highest receives 1; lowest receives length(C); tied values receive same rank)
        nn = observation_year;
        P = (rr)./(nn+1);
        RP = 1./P;
        RP = sort(RP,'ascend');
        MRI = sort(maximas,'ascend');
    end

end