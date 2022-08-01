function [r_mean,r_var,pd_mean,pd_var,pd_auto,dc_var,dd_var] = msm_calculate_moments(allvariables_levels)
%Calculates moments for msm estimation



%% Equity return
%Get equity returns
[r] = get_returns(allvariables_levels(:,1),allvariables_levels(:,4),allvariables_levels(:,5));

%Calculate mean equity return
r_mean = mean(r);
r_var = var(r,1);

%% Price Dividend Ratio
%Calculate Price-Dividend ratio
dimensions_allvariables_levels = size(allvariables_levels);
n = dimensions_allvariables_levels(1);
pd = zeros(n,1);

for i = 4:n
    pd(i) = allvariables_levels(i,1)/(allvariables_levels(i,4)+allvariables_levels(i-1,4)+allvariables_levels(i-2,4)+allvariables_levels(i-3,4));

end

%Calculate moments
pd = pd(4:end); 
pd_mean = mean(pd);
pd_std = std(pd,1);
pd_var = pd_std^2;
pd_cov = cov(pd(2:end),pd(1:end-1),1);
pd_auto = pd_cov(2,1);

%% Macro Variables
%Calculate change in macro variables
dc = zeros(n,1);
dd = zeros(n,1);

for i = 2:n
    dc(i) = log(allvariables_levels(i,3)) - log(allvariables_levels(i-1,3)); %Change in consumption
    dd(i) = log(allvariables_levels(i,4)) - log(allvariables_levels(i-1,4)); %Change in dividends

end


%Calculate moments

%Consumption
dc = dc(2:n);
dc_var = var(dc,1);

%Dividends
dd = dd(2:n);
dd_var = var(dd,1);


end