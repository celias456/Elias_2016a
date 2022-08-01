function [data] = statistics_macro(allvariables_levels)
%Calculates statistics on macro data (log growth rates)

%Returns a 1x10 vector
%   dc_std: change in consumption std
%   dc_var: change in consumption variance
%   c_corr: change in consumption first order autocorrelation
%   dd_std: change in dividends std
%   dd_var: change in dividends variance
%   d_corr: change in dividends first order autocorrelation
%   dy_std: change in output std
%   y_corr: change in output first order autocorrelation
%   di_std: change in investment std
%   i_corr: investment first order autocorrelation
% Funtions used by this code
% (1) RETURNS.m 


%Calculate change in macro variables
dimensions_allvariables_levels = size(allvariables_levels);
n = dimensions_allvariables_levels(1);
dc = zeros(n,1);
dd = zeros(n,1);
dy = zeros(n,1);
dk = zeros(n,1);
di = zeros(n,1);
dgp = zeros(n,1);

for i = 2:n
    dc(i) = log(allvariables_levels(i,3)) - log(allvariables_levels(i-1,3)); %Change in consumption
    dd(i) = log(allvariables_levels(i,4)) - log(allvariables_levels(i-1,4)); %Change in dividends
    
    dy(i) = log(allvariables_levels(i,2)) - log(allvariables_levels(i-1,2)); %Change in output
    dk(i) = log(allvariables_levels(i,6)) - log(allvariables_levels(i-1,6)); %Change in capital
    di(i) = log(allvariables_levels(i,8)) - log(allvariables_levels(i-1,8)); %Change in investment
    dgp(i) = log(allvariables_levels(i,9)) - log(allvariables_levels(i-1,9)); %Change in gross profit

end


%Calculate moments and first-order auto-correlations

%Consumption
dc = dc(2:n);
dc_std = 200*std(dc,1);
dc_var = var(dc,1);
c_corr_matrix = corr(dc(2:end),dc(1:end-1));
c_corr = c_corr_matrix(1,1);

%Dividends
dd = dd(2:n);
dd_std = 200*std(dd,1);
dd_var = var(dd,1);
d_corr_matrix = corr(dd(2:end),dd(1:end-1));
d_corr = d_corr_matrix(1,1);

%Output
dy = dy(2:n);
dy_std = 200*std(dy,1);
y_corr_matrix = corr(dy(2:end),dy(1:end-1));
y_corr = y_corr_matrix(1,1);

% %Capital
% dk = dk(2:n); 
% dk_std = 200*std(dk);
% k_corr_matrix = corr(dk(2:end),dk(1:end-1));
% k_corr = k_corr_matrix(1,1);

%Investment
di = di(2:n); 
di_std = 200*std(di,1);
i_corr_matrix = corr(di(2:end),di(1:end-1));
i_corr = i_corr_matrix(1,1);

% %Gross Profit
% dgp = dgp(2:n);
% dgp_std = 200*std(dgp);
% gp_corr_matrix = corr(dgp(2:end),dgp(1:end-1));
% gp_corr = gp_corr_matrix(1,1);


data = [dc_std,dc_var,c_corr,dd_std,dd_var,d_corr,dy_std,y_corr,di_std,i_corr];