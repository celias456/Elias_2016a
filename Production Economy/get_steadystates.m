function [ss] = get_steadystates(alpha,beta,delta,psi)
%Calculates value of steady states given parameters of the model
%Returned matrix is steady state value of equity price, output,
%consumption, dividends, bond price, capital, productivity shock,
%investment, and gross profit

k = psi^(-1/(1-alpha)); %Capital
y = k^alpha; %Output
c = k^alpha - delta*k; %Consumption
i = delta*k; %Investment
d = alpha*k^alpha - delta*k; %Dividends
pb = beta; %Bond Price
x = alpha*y; %Gross Profits
p = k; %Equity Price
z = 1; %Productivity Shock

ss = [p;y;c;d;pb;k;z;i;x];