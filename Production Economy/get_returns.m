function [R,Rf,ep,R_gross,Rf_gross,count] = get_returns(p,d,pb)
%Calculate equity and risk free percent and gross returns

T = length(p);
count = 0;
R = zeros(T,1);
Rf = zeros(T,1);
ep = zeros(T,1);

for t = 2:T
    tm1 = t-1;
    R(t) = (((p(t) + d(t))/p(tm1))-1)*100;
    
    if pb(tm1) == 0
       pb(tm1) = 1e-323;
       count = count + 1;
    end
    
    Rf(t) = ((1/pb(tm1))-1)*100;
    ep(t) = (R(t) - Rf(t));
end

R = R(2:end); %equity returns percent
Rf = Rf(2:end); %risk free rates percent
ep = ep(2:end); %equity premium percent

%Calculate gross returns (1+r)
R_gross = (R/100)+1;
Rf_gross = (Rf/100)+1;