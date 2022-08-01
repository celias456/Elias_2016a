function [B] = numderiv(estimates,tol,sT,EM)
%Calculates Numrical derivatives of the simulated moments
%Numerical derivatives are calculated using a two-point estimation method

mu = estimates(1);
eshock_std = estimates(2);
eshock_ro = estimates(3);
gain = estimates(4);

B = zeros(4,4);

all = g_st([mu,eshock_std,eshock_ro,gain],sT,EM);

%B is a 4x4 matrix where the row represents the moment and the column
%represents the parameter

%Calculate numerical derivatives for mu
mu_ub_all = g_st([mu+tol,eshock_std,eshock_ro,gain],sT,EM);
for count = 1:4
    B(count,1) = (mu_ub_all(count)-all(count))/(tol);
end

%Calculate numerical derivatives for eshock_std
eshock_std_ub_all = g_st([mu,eshock_std+tol,eshock_ro,gain],sT,EM);
for count = 1:4
    B(count,2) = (eshock_std_ub_all(count)-all(count))/(tol);
end

%Calculate numerical derivatives for eshock_ro
eshock_ro_ub_all = g_st([mu,eshock_std,eshock_ro+tol,gain],sT,EM);
for count = 1:4
    B(count,3) = (eshock_ro_ub_all(count)-all(count))/(tol);
end

%Calculate numerical derivatives for gain
gain_ub_all = g_st([mu,eshock_std,eshock_ro,gain+tol],sT,EM);
for count = 1:4
    B(count,4) = (gain_ub_all(count)-all(count))/(tol);
end

end
