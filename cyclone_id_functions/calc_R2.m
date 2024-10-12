function R2 = calc_R2(x,dist)
n=length(x);
x = sort(x,"ascend");
p = ((1:n) - 0.5) ./ n;
q_dist = icdf(dist,p)';
R2 = 1 - sum((x - q_dist).^2)/sum((x - mean(x)).^2);
end