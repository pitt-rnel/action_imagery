function[lower_bound, upper_bound, sort_boot, rand_boot] = boot_bounds(nboot,bootfun,data,lower_perc,upper_perc)
% Helper function to perform bootstrapping 

if lower_perc > 1 || upper_perc > 1
    lower_perc = lower_perc./100;
    upper_perc = upper_perc./100;
end

rand_boot = bootstrp(nboot,bootfun,data);
sort_boot = sortrows(rand_boot);

lower_upper = prctile(rand_boot,100*[lower_perc upper_perc],1);

lower_bound = lower_upper(1,:)';
upper_bound = lower_upper(2,:)';

end

