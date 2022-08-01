function [lom,pf,count,Rtest1,Rtest2] = lom_ree(a,b,ro,gamma,Pbar,Ztm2,shock,mu,T,eshock,eshock_wn) 
% Returns Law of motion for all variables under a REE
% Returns a Tx7 matrix lom containing the law of motion under a HEE for
% variables
%    lom(:,1) = Equity price
%    lom(:,2) = Output
%    lom(:,3) = Consumption
%    lom(:,4) = Dividends
%    lom(:,5) = Bond Price
%    lom(:,6) = Capital
%    lom(:,7) = shock
%    lom(:,8) = investment
%    lom(:,9) = gross profits
 

%Initialize storage 
lom = zeros(T+2,9); %Stores law of motion of all variables
expshocks = zeros(T+2,2); %Stores law of motion of expectation shocks

phibar_c1 = Pbar(1);
phibar_d = Pbar(2);
phibar_c2 = Pbar(1);

%Set initial values
Ztm1 = Ztm2*ro + shock(2); %Initial value of shock
Dtm1 = Ztm1;
lom(2,7) = Ztm1; %Store initial value of shock
lom(2,4) = Dtm1; %Store initial value of dividends
expshocks(2,1) = eshock(1,1)*(eshock_wn(1,1)*eshock(1,4) + eshock_wn(2,1)); %Store initial value of expectation shock for agent 1
expshocks(2,2) = eshock(2,1)*(eshock_wn(1,2)*eshock(2,4) + eshock_wn(2,2)); %Store initial value of expectation shock for agent 2
EP_0 = phibar_d*ro*Dtm1 + expshocks(2,1); %Agent expectations 
lom(2,1) = a*EP_0 + b*Dtm1; %Store initial value of price


for t = 3:(T+2)
    
   tm1 = t-1;
   Z_tm1 = lom(tm1,7); %Productivity shock in previous period
   Z_t = Z_tm1*ro + shock(t); %Current value of productivity shock
   D_t = Z_t; %Current value of dividends
   expshocks(t,1) = eshock(1,1)*(expshocks(tm1,1)*eshock(1,4) + eshock_wn(t,1)); %Current value of expectation shock for agent 1
   expshocks(t,2) = eshock(2,1)*(expshocks(tm1,2)*eshock(2,4) + eshock_wn(t,2)); %Current value of expectation shock for agent 2
   
   %Obtain law of motion for variables under REE
   EP_tp1 = phibar_d*ro*D_t + expshocks(t,1); %Agent expectations
   P_t = a*EP_tp1 + b*D_t; 
   
   %Obtain law of motion for other endogeneous variables in the model
   lom(t,1) = P_t;                                                 %Price
   lom(t,2) = gamma(1,1)*P_t + gamma(1,2)*D_t + gamma(1,3)*EP_tp1; %Output
   lom(t,3) = gamma(2,1)*P_t + gamma(2,2)*D_t + gamma(2,3)*EP_tp1; %Consumption
   lom(t,4) = D_t;                                                 %Dividends
   lom(t,5) = gamma(4,1)*P_t + gamma(4,2)*D_t + gamma(4,3)*EP_tp1; %Bond Price
   lom(t,6) = gamma(5,1)*P_t + gamma(5,2)*D_t + gamma(5,3)*EP_tp1; %Capital
   lom(t,7) = Z_t;                                                 %Shock
   lom(t,8) = gamma(6,1)*P_t + gamma(6,2)*D_t + gamma(6,3)*EP_tp1; %Investment
   lom(t,9) = gamma(7,1)*P_t + gamma(7,2)*D_t + gamma(7,3)*EP_tp1; %Gross Profit
   
   %Obtain law of motion for agent forecasts
   lom(t,10) = 0; %Previous period price
   lom(t,11) = 0; %Agent 1 previous period prediction of price
   lom(t,12) = 0; %Agent 2 previous period prediction of price
      
end

lom = lom(3:T+2,:);
pf = 0;
count = 0;
Rtest1 = 0;
Rtest2 = 0;

end