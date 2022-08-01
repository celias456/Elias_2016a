function [level] = get_levels(ss,T,allvariables,burn)
%Returns levels of variables 
%Takes as input series of log deviations from steady state and returns
%levels 

begin = burn+1;
keep = T-burn;
level = zeros(keep,12);

level(:,1) = ss(1)*exp(allvariables(begin:T,1)); %Price
level(:,2) = ss(2)*exp(allvariables(begin:T,2)); %Output
level(:,3) = ss(3)*exp(allvariables(begin:T,3)); %Consumption
level(:,4) = ss(4)*exp(allvariables(begin:T,4)); %dividends
level(:,5) = ss(5)*exp(allvariables(begin:T,5)); %Bond Price
level(:,6) = ss(6)*exp(allvariables(begin:T,6)); %Capital
level(:,7) = ss(7)*exp(allvariables(begin:T,7)); %Shock
level(:,8) = ss(8)*exp(allvariables(begin:T,8)); %Investment
level(:,9) = ss(9)*exp(allvariables(begin:T,9)); %Gross Profit

%Prediction error
level(:,10) = ss(1)*exp(allvariables(begin:T,10)); %previous period price
level(:,11) = ss(1)*exp(allvariables(begin:T,11)); %Agent 1 previous period prediction
level(:,12) = ss(1)*exp(allvariables(begin:T,12)); %Agent 2 previous period prediction
