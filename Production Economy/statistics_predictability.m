function [data] = statistics_predictability(r_gross,rf_gross,pdratio)
%Calculates predictability statistics

%Returns a 1x6 vector 'data'
%data(1): slope coefficient for 1 year ahead regression
%data(2): = 1 if coefficient is significant
%data(3): adjusted R-squared for 1 year ahead regression
%data(4): slope coefficient for 2 year ahead regression
%data(5): = 1 if coefficient is significant
%data(6): adjusted R-squared for 2 year ahead regression
%data(7): slope coefficient for 4 year ahead regression
%data(8): = 1 if coefficient is significant
%data(9): adjusted R-squared for 4 year ahead regression

%Find 1,2,and 4 year ahead cumulative returns
[r_cum_1,r_cum_2,r_cum_4] = get_returns_cumulative(r_gross); %equity 
[rf_cum_1,rf_cum_2,rf_cum_4] = get_returns_cumulative(rf_gross); %risk free rate

%Calculate 1,2, and 4 year ahead equity premium in natural units
ep_1 = (r_cum_1 - rf_cum_1)/ 100; 
ep_2 = (r_cum_2 - rf_cum_2)/ 100; 
ep_4 = (r_cum_4 - rf_cum_4)/ 100; 

%Take natural log of price dividend ratio
pdratio = log(pdratio);

%Find sample standard deviation of log price dividend ratio
pdratio_sd = std(pdratio,1);

%Transform log price dividend ratio into standard deviation units
pdratio = pdratio/pdratio_sd;

%Regress the 1 yr ahead equity premium on the transformed price
%dividend ratio
pdratio_1 = pdratio(1:end-3);
ep_1 = ep_1(3:end);
t_1 = length(pdratio_1);
l_1 = floor(4*(t_1/100)^(2/9));
c_1 = ones(t_1,1);
x_1 = [c_1,pdratio_1];
reg_1 = nwest(ep_1,x_1,l_1);

beta_1 = reg_1.beta(2);
if reg_1.tstat(2) <= -1.96
    beta_t_1 = 1;
else
    beta_t_1 = 0;
end
rsqr_1 = reg_1.rbar;

% %Regress the 2 yr ahead equity premium on the transformed price
% %dividend ratio
pdratio_2 = pdratio(1:end-7);
ep_2 = ep_2(3:end);
t_2 = length(pdratio_2);
l_2 = floor(4*(t_2/100)^(2/9));
c_2 = ones(t_2,1);
x_2 = [c_2,pdratio_2];
reg_2 = nwest(ep_2,x_2,l_2);

beta_2 = reg_2.beta(2);
if reg_2.tstat(2) <= -1.96
    beta_t_2 = 1;
else
    beta_t_2 = 0;
end
rsqr_2 = reg_2.rbar;

% %Regress the 4 yr ahead equity premium on the transformed price
% %dividend ratio
pdratio_4 = pdratio(1:end-15);
ep_4 = ep_4(3:end);
t_4 = length(pdratio_4);
l_4 = floor(4*(t_4/100)^(2/9));
c_4 = ones(t_4,1);
x_4 = [c_4,pdratio_4];
reg_4 = nwest(ep_4,x_4,l_4);

beta_4 = reg_4.beta(2);
if reg_4.tstat(2) <= -1.96
    beta_t_4 = 1;
else
    beta_t_4 = 0;
end
rsqr_4 = reg_4.rbar;

data = [beta_1,beta_t_1,rsqr_1,beta_2,beta_t_2,rsqr_2,beta_4,beta_t_4,rsqr_4];

end

