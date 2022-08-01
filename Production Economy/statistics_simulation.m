function [stats_returns,stats_pdratio,stats_predictability,stats_macro,stats_errors,stats_moments,stats_prediction_errors] = statistics_simulation(data_returns,data_pdratio,data_predictability,data_macro,data_errors,data_predict_errors)
%Calculates statistics from a number of simulations


%Calculate return medians
returns = median(data_returns,1);

%Calculate price-dividend ratio medians
pdratio = median(data_pdratio,1);

%Predictability
predict = median(data_predictability,1); %Calculate medians
dim = size(data_predictability);
n = dim(1,1); 
totals = sum(data_predictability,1);
t_1 = (totals(2)/n)*100;
t_2 = (totals(5)/n)*100;
t_4 = (totals(8)/n)*100;

%Calculate macro data medians
macro = median(data_macro,1);

%Calculate prediction error medians;
prediction_errors = median(data_predict_errors,1);

pf1_totaltimes = sum(data_errors(:,1));
pf2_totaltimes = sum(data_errors(:,2));
pf1_percent = (sum(data_errors(:,3)))/n;
pf2_percent = (sum(data_errors(:,4)))/n;
r1 = sum(data_errors(:,5));
r2 = sum(data_errors(:,6));

stats_returns = struct('ER_mean',returns(1),'ER_std',returns(2),'RF_mean',returns(3),'RF_std',returns(4),'EP_mean',returns(5),'EP_std',returns(6));
stats_pdratio = struct('PD_mean',pdratio(1),'PD_std',pdratio(2),'PD_corr',pdratio(4),'Div_std',macro(4),'Div_corr',macro(6));
stats_predictability = struct('beta_1',predict(1),'t_1',t_1,'R_1',predict(3),'beta_2',predict(4),'t_2',t_2,'R_2',predict(6),'beta_4',predict(7),'t_4',t_4,'R_4',predict(9));
stats_macro = struct('C_std',macro(1),'C_corr',macro(3),'Y_std',macro(7),'Y_corr',macro(8),'I_std',macro(9),'I_corr',macro(10));
stats_errors = struct('pf1_totaltimes',pf1_totaltimes,'pf2_totaltimes',pf2_totaltimes,'pf1_percent',pf1_percent,'pf2_percent',pf2_percent,'r1',r1,'r2',r2);
stats_moments = struct('ER_mean',returns(1),'ER_var',returns(7),'PD_mean',pdratio(1),'PD_var',pdratio(3),'PD_auto',pdratio(5),'C_var',macro(2),'D_var',macro(5));
stats_prediction_errors = struct('Agent_1', prediction_errors(1), 'Agent_2', prediction_errors(2));
end