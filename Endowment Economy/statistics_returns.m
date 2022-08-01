function [data,r_gross,rf_gross] = statistics_returns(allvariables_levels)
%Calculates return statistics

%Returns a 1x8 vector
%   r_mean = equity return mean
%   r_std = equity return sd
%   rf_mean = risk free rate mean
%   rf_std = risk free rate sd
%   ep_mean = equity premium mean
%   ep_std = equity premium sd
%   r_gross = equity gross return
%   rf_rf gross return


%Get returns and equity premium
[r,rf,ep,r_gross,rf_gross] = get_returns(allvariables_levels(:,1),allvariables_levels(:,4),allvariables_levels(:,5));

%Calculate average rate of return
r_mean = mean(r);
rf_mean = mean(rf);
ep_mean = mean(ep);

%Calculate standard deviation of returns
r_std = std(r,1);
r_var = r_std^2;
rf_std = std(rf,1);
ep_std = std(ep,1);


data = [r_mean,r_std,rf_mean,rf_std,ep_mean,ep_std,r_var];

end