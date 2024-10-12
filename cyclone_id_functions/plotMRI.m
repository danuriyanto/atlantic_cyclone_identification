function [RP,MRI] =  plotMRI(maximas,observation_year)
    [C,~,ic] = unique(maximas,'sorted'); % ic are ranks from lowest to highest ; C are unique values
    rr=(1+max(ic)-ic);  % r: rank (highest receives 1; lowest receives length(C); tied values receive same rank)
    nn = observation_year;
    P = (rr)./(nn+1);
    RP = 1./P;
    RP = sort(RP,'ascend');
    MRI = sort(maximas,"ascend");
end