function [dist] = GEV_LS_or_MLE(extremes_vector)

[dist_LS] = LS_GEV_limit(extremes_vector);
[dist_MLE] = MLE_GEV_limit(extremes_vector);

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