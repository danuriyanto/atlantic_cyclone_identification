function [dist] = GPD_LS_or_MLE(extremes_vector,thres)

[dist_LS] = LS_GPD_limit(extremes_vector,thres);
[dist_MLE] = MLE_gpd_limit(extremes_vector,thres);

R2_LS = calc_R2(extremes_vector,dist_LS);
R2_MLE = calc_R2(extremes_vector,dist_MLE);

if R2_MLE > R2_LS
    sprintf('R2 MLE = %0.03f R2 LS = %0.03f, we chose MLE',R2_MLE,R2_LS)
    dist = dist_MLE;
elseif R2_MLE == R2_LS
    dist = dist_LS;
    sprintf('R2 MLE = %0.03f R2 LS = %0.03f, all similar',R2_MLE,R2_LS)

else
    dist = dist_LS;
    sprintf('R2 MLE = %0.03f R2 LS = %0.03f, we chose LS',R2_MLE,R2_LS)
end


end