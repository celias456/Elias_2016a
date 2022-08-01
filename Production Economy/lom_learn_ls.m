function [lom,pf1,pf2,count1,count2,Rtest1,Rtest2] = lom_learn_ls(a1,a2,b1,ro,gamma,P,Mbar1,Mbar2,Ptm2,Ztm2,shock,mu,T,eshock,eshock_wn,constant_2) 
% Returns Law of motion for all variables under RLS AL
% Returns a Tx7 matrix lom containing the law of motion under RLS learning for
% variables
%    low(:,1) = Equity price
%    low(:,2) = Output
%    low(:,3) = Consumption
%    low(:,4) = Dividends
%    low(:,5) = Bond Price
%    low(:,6) = Capital
%    low(:,7) = Productivity shock
%    lom(:,8) = Investment
%    lom(:,9) = Gross profits

minmu = 1 - mu;

%Create storage matrices
lom = zeros(T+2,12); %Stores law of motion of all variables
Phi = zeros(3,T+1); %Store values for the parameter updates (agent 1) 
Theta = zeros(2,T+1); %Store values for the parameter updates (agent 2) 
R1 = zeros(3,3,T+1); %Store values for R matrix (agent 1)

%Store values for R matrix (agent 2)
if constant_2 == 1
    R2 = zeros(2,2,T+1);
else
    R2 = zeros(1,T+1);
end

expshocks = zeros(T+2,2); %Stores law of motion of expectation shocks

%Initialize variables
R1(:,:,1) = Mbar1; %Initial value for R matrix (agent 1)

%Initial value for R matrix (agent 2)
if constant_2 == 1
    R2(:,:,1) = Mbar2;
else
    R2(1) = Mbar2; 
end

Phi(:,1) = P(1:3); %Initial values for parameters which are updated each period (agent 1)
Theta(:,1) = P(4:5); %Initial values for parameters which are updated each period (agent 2)
PhiC1_0 = Phi(1,1); %Initial parameter value for constant (agent 1)
PhiP_0 = Phi(2,1); %Initial parameter value on lagged price (agent 1)
PhiZ_0 = Phi(3,1); %Initial parameter value on current shock (agent 1)
PhiC2_0 = Theta(1,1); %Initial parameter value for constant (agent 2)
theta_0 = Theta(2,1); %Initial parameter value on current shock (agent 2)
Ztm1 = Ztm2*ro + shock(2); %Initial value of shock
expshocks(2,1) = eshock(1,1)*(eshock_wn(1,1)*eshock(1,4) + eshock_wn(2,1)); %Store initial value of expectation shock for agent 1
expshocks(2,2) = eshock(2,1)*(eshock_wn(1,2)*eshock(2,4) + eshock_wn(2,2)); %Store initial value of expectation shock for agent 2

%Initial expectations of next period equity price of agents
EP1_0 = PhiC1_0 + PhiP_0*PhiC1_0 + (PhiP_0^2)*Ptm2 + (PhiP_0*PhiZ_0 + PhiZ_0*ro)*Ztm1 + expshocks(2,1); %Agent 1
EP2_0 = PhiC2_0*(1+theta_0) + (theta_0^2)*Ptm2 + expshocks(2,2);%Agent 2
EP_0 = mu*EP1_0 + (minmu)*EP2_0; %Average expectations in economy

%Store initial values for equity price and shock
lom(1,1) = Ptm2; %Store initial value of price
lom(2,1) = a1*EP_0 + a2*Ptm2 + b1*Ztm1; %Store second period equity price (From ALM)
lom(2,7) = Ztm1; %Store initial value of shock

%Initialize counters 
pf1 = 0; %Counts number of times projection facility is used (agent 1).
pf2 = 0; %Counts number of times projection facility is used (agent 2).
Rtest1 = 0; %Counts number of times the R1 matrix is not invertible
Rtest2 = 0; %Counts number of times the R2 matrix is not invertible

%Begin RLS learning recursion

for t = 3:(T+2)
       
   %Set values for recursion
   tm2 = t-2; 
   tm1 = t-1;
   P_tm2 = lom(tm2,1); %Equity price two periods previously
   P_tm1 = lom(tm1,1); %Equity price one period previously
   Z_tm1 = lom(tm1,7); %Shock in previous period
   gain = 1/tm2;
   Phi_tm1 = Phi(:,tm2); %Values of agent 1 coefficient estimates one period previously
   PhiC1_tm1 = Phi_tm1(1); %Agent 1 constant one period previously
   PhiP_tm1 = Phi_tm1(2); %Agent 1 coefficient on lagged equity price one period previously
   PhiZ_tm1 = Phi_tm1(3); %Agent 1 coefficient on productivity shock one period previously
   
   Theta_tm1 = Theta(:,tm2); %Values of agent 2 coefficient estimates one period previously
   PhiC2_tm1 = Theta_tm1(1); %Agent 2 constant one period previously
   theta_tm1 = Theta_tm1(2); %Agent 2 coefficient on lagged equity price one period previously
   
   R1_tm1 = R1(:,:,tm2); %Value of agent 1 R matrix one period previously
   
   %Value of agent 2 R matrix one period previously   
   if constant_2 == 1
       R2_tm1 = R2(:,:,tm2);
   else
       R2_tm1 = R2(tm2); 
   end

   Phat1_tm1 = PhiC1_tm1 + PhiP_tm1*P_tm2 + PhiZ_tm1*Z_tm1; %Agent 1 most recent predicted value for equity price
   prederror1 = P_tm1 - Phat1_tm1; %Agent 1 most recent prediction error
   
   Phat2_tm1 = PhiC2_tm1 + theta_tm1*P_tm2; %Agent 2 most recent predicted value for equity price
   prederror2 = P_tm1 - Phat2_tm1; %Agent 2 most recent prediction error

   
   %Agents update parameters by adding P_tm1 to the information set and run
   %a LS regression of P_tm1 on a constant, P_tm2, and X_tm1.

   %Update the R matrix - Agent 1
   R1_t = R1_tm1 + gain * ( ([1;P_tm2;Z_tm1] * [1,P_tm2,Z_tm1]) - R1_tm1); 
   R1inv = (R1_t)^-1;
   %Update the coefficient matrix - Agent 1
   Phi_t = Phi_tm1 + gain * R1inv * ([1;P_tm2;Z_tm1]*(prederror1)); 

   %Update the R matrix - Agent 2
   if constant_2 == 1
       R2_t = R2_tm1 + gain * ( ([1;P_tm2] * [1,P_tm2]) - R2_tm1);
       R2inv = (R2_t)^-1;
       %Update the coefficient matrix - Agent 2
       Theta_t = Theta_tm1 + gain * R2inv * ([1;P_tm2]*(prederror2));
   else
       R2_t = R2_tm1 + gain * ( (P_tm2 * P_tm2) - R2_tm1); 
       R2inv = (R2_t)^-1;
       %Update the coefficient matrix - Agent 2
       Theta_t(2,1) = theta_tm1 + gain * R2inv * (P_tm2*(prederror2));
   end
   
   PhiC1_t = Phi_t(1,1);
   PhiP_t = Phi_t(2,1);
   PhiZ_t = Phi_t(3,1);
   
   PhiC2_t = Theta_t(1,1);
   theta_t = Theta_t(2,1);

   %This is a check to make sure that R1 and R2 are invertible
   %Agent 1
   if R1inv(1,1) == inf || R1inv(1,2) == inf || R1inv(1,3) == inf || R1inv(2,1) == inf || R1inv(2,2) == inf || R1inv(2,3) == inf || R1inv(3,1) == inf || R1inv(3,2) == inf || R1inv(3,3) == inf
      PhiC1_t = PhiC1_tm1; 
      PhiP_t = PhiP_tm1;
      PhiZ_t = PhiZ_tm1;
      Rtest1 = Rtest1 + 1;
   end
   
   %Agent 2
   if constant_2 == 1
       if R2inv(1,1) == inf || R2inv(1,2) == inf || R2inv(2,1) == inf || R2inv(2,2) == inf
           PhiC2_t = PhiC2_tm1;
           theta_t = theta_tm1;
           Rtest2 = Rtest2 + 1;
       end
       
   else
       if R2inv == inf
           PhiC2_t = PhiC2_tm1;
           theta_t = theta_tm1;
           Rtest2 = Rtest2 + 1;
       end
   end

   %Projection Facility
   %Agent 1
   if abs(PhiP_t) >= 1
      PhiC1_t = PhiC1_tm1;
      PhiP_t = PhiP_tm1;
      PhiZ_t = PhiZ_tm1;
      pf1 = pf1 + 1;
   end
   %Agent 2
   if abs(theta_t) >= 1
      PhiC2_t = PhiC2_tm1;
      theta_t = theta_tm1;
      pf2 = pf2 + 1;
   end
   
   %Update Phi, R1, Theta, and R2 with values of new parameters
   Phi(1,tm1) = PhiC1_t;
   Phi(2,tm1) = PhiP_t;
   Phi(3,tm1) = PhiZ_t;
   R1(:,:,tm1) = R1_t;
   
   Theta(1,tm1) = PhiC2_t;
   Theta(2,tm1) = theta_t;
   if constant_2 == 1
       R2(:,:,tm1) = R2_t;
   else
       R2(tm1) = R2_t;
   end
   
   %The value of the current period shock (Z_t) is realized
   Z_t = Z_tm1*ro + shock(t);
   
   %Expectation shocks are realized
   expshocks(t,1) = eshock(1,1)*(expshocks(tm1,1)*eshock(1,4) + eshock_wn(t,1)); %Current value of expectation shock for agent 1
   expshocks(t,2) = eshock(2,1)*(expshocks(tm1,2)*eshock(2,4) + eshock_wn(t,2)); %Current value of expectation shock for agent 2
   
   %Agents form expectations of the future price of equity 
   EP1_tp1 = PhiC1_t + PhiP_t*PhiC1_t + (PhiP_t^2)*P_tm1 + (PhiP_t*PhiZ_t+PhiZ_t*ro)*Z_t + expshocks(t,1); %Agent 1
   EP2_tp1 = PhiC2_t*(1+theta_t) + (theta_t^2)*P_tm1 + expshocks(t,2); %Agent 2
   EP_tp1 = mu*EP1_tp1 + (minmu)*EP2_tp1; %Average expectations in economy
   
   %P_T is realized from the reduced form of the model (or from ALM)
   P_t = a1*EP_tp1 + a2*P_tm1 + b1*Z_t;
   
   
   %Obtain law of motion for dividends and bond price under adaptive learning
   lom(t,1) = P_t; %Price
   lom(t,2) = gamma(1,1)*P_t + gamma(1,2)*P_tm1 + gamma(1,3)*Z_t + gamma(1,4)*EP_tp1; %Output
   lom(t,3) = gamma(2,1)*P_t + gamma(2,2)*P_tm1 + gamma(2,3)*Z_t + gamma(2,4)*EP_tp1; %Consumption
   lom(t,4) = gamma(3,1)*P_t + gamma(3,2)*P_tm1 + gamma(3,3)*Z_t + gamma(3,4)*EP_tp1; %Dividends
   lom(t,5) = gamma(4,1)*P_t + gamma(4,2)*P_tm1 + gamma(4,3)*Z_t + gamma(4,4)*EP_tp1; %Bond Price
   lom(t,6) = gamma(5,1)*P_t + gamma(5,2)*P_tm1 + gamma(5,3)*Z_t + gamma(5,4)*EP_tp1; %Capital
   lom(t,7) = Z_t; %Shock
   lom(t,8) = gamma(6,1)*P_t + gamma(6,2)*P_tm1 + gamma(6,3)*Z_t + gamma(6,4)*EP_tp1; %Investment
   lom(t,9) = gamma(7,1)*P_t + gamma(7,2)*P_tm1 + gamma(7,3)*Z_t + gamma(7,4)*EP_tp1; %Gross Profit
   
   %Obtain law of motion for agent forecasts
   lom(t,10) = P_tm1; %Previous period price
   lom(t,11) = Phat1_tm1; %Agent 1 previous period prediction of price
   lom(t,12) = Phat2_tm1; %Agent 2 previous period prediction of price
   
      
end
   
%Drop the first two values in all the series
lom = lom(3:T+2,:);

if pf1 > 0
   count1 = 1;
else
   count1 = 0;
end

if pf2 > 0
   count2 = 1;
else
   count2 = 0;
end

end

   
    

