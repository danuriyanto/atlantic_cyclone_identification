function  bootstrap_result = bootstrapping(extremes_vector,distribution,thres,num_bootstrap)
boot_k = zeros(num_bootstrap,1);
boot_sig = zeros(num_bootstrap,1);
boot_thres = zeros(num_bootstrap,1);

for i = 1:num_bootstrap
    % Generate a bootstrap sample
    bootstrap_sample = random(distribution,size(extremes_vector));
    
    % Fit GPD to the bootstrap sample
    bootstrap_params = LS_GPD_limit2(bootstrap_sample,thres);

    % store the parameters
    boot_k(i)          = bootstrap_params.k;
    boot_sig(i)        = bootstrap_params.sigma;
    boot_thres(i)      = bootstrap_params.theta;
end

bootstats_k = mean(boot_k);
bootstats_sig = mean(boot_sig);
bootstats_thres = mean(boot_thres);

bootstrap_result = struct();
bootstrap_result.description        = {'k - shape' 'sigma - scale' 'theta - threshold/location'};
bootstrap_result.bootstrap_params   = [bootstats_k bootstats_sig bootstats_thres];
bootstrap_result.fitting_params     = [distribution.k distribution.sigma distribution.theta];
bootstrap_result.standard_error     = [std(boot_k) std(boot_sig) std(boot_thres)];
bootstrap_result.param_ci_95prct    = [prctile(boot_k, [2.5 97.5]); prctile(boot_sig, [2.5 97.5]); prctile(boot_thres, [2.5 97.5])];

end
