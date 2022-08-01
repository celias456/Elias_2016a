%% Equity Price Graph
load equityprice.csv
equityprice_data = equityprice;

x1 = equityprice_data(:,1);
y1 = equityprice_data(:,2);
y2 = equityprice_data(:,3);

% plot(x1,y1,'k',x1,y2,'--k')
% xlabel('Period')
% ylabel('Equity Price')


figure
subplot(2,1,1)
plot(x1,y1)
xlabel('Period')
ylabel('Equity Price')
title('RE1')

subplot(2,1,2)
plot(x1,y2)
xlabel('Period')
ylabel('Equity Price')
title('ALE')


%% P-D Ratio Graph

load pdratio.csv
pdratio_data = pdratio;

x1 = pdratio_data(:,1);
y1 = pdratio_data(:,2);
y2 = pdratio_data(:,3);

plot(x1,y1,'k',x1,y2,'--k')
xlabel('Period')
ylabel('P-D Ratio')

legend('RE1','ALE')

%% Forecast Graph

load forecast.csv
forecast_data = forecast;

x1 = forecast_data(:,1);
y1 = forecast_data(:,2);
y2 = forecast_data(:,3);
y3 = forecast_data(:,4);

plot(x1,y1,'k',x1,y2,'--k',x1,y3,':k')
xlabel('Period')
% ylabel('Equity Price')

legend('Active Trader Forecast','Passive Investor Forecast','Equity Price')

%% Beliefs

load beliefs.csv
beliefs_data = beliefs;

x1 = beliefs_data(:,1);
y1 = beliefs_data(:,2);
y2 = beliefs_data(:,3);
y3 = beliefs_data(:,4);


figure
subplot(3,1,1)
plot(x1,y1)
xlabel('Period')
ylabel('Constant')

subplot(3,1,2)
plot(x1,y2)
xlabel('Period')
ylabel('Lagged Price')

subplot(3,1,3)
plot(x1,y3)
xlabel('Period')
ylabel('Productivity Shock')

