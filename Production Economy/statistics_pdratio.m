function [data,pd] = statistics_pdratio(allvariables_levels)
%Calculates price dividend ratio statistics

%Returns a 1x4 vector
%   pd_mean: price-dividend ratio mean
%   pd_std: price-dividend ratio std
%   pd_var: price-dividend variance
%   pd_corr: price dividend first order autocorrelation

%Returns a (n-3) x 1 vector of the price dividend ratio


%Calculate Price-Dividend ratio and change in macro variables
dimensions_allvariables_levels = size(allvariables_levels);
n = dimensions_allvariables_levels(1);
pd = zeros(n,1);


for i = 4:n
    pd(i) = allvariables_levels(i,1)/(allvariables_levels(i,4)+allvariables_levels(i-1,4)+allvariables_levels(i-2,4)+allvariables_levels(i-3,4));

end

%Calculate moments and first-order auto-correlations
pd = pd(4:end);
k = length(pd); 
pd_mean = mean(pd);
pd_std = std(pd,1);
pd_var = pd_std^2;
cc = corr(pd(2:k),pd(1:k-1));
pd_corr = cc(1,1);
pd_cov = cov(pd(2:end),pd(1:end-1),1);
pd_auto = pd_cov(2,1);



data = [pd_mean,pd_std,pd_var,pd_corr,pd_auto];
