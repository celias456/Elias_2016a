function [B] = msm_numderiv(seedstate,estimates,tol,sT,EM,expectationshock,constant_2,n)
%Calculates Numrical derivatives of the simulated moments
%Numerical derivatives are calculated using a two-point estimation method

%Returns a 5x5 matrix 'B' where the row represents the moment and the column
%represents the parameter

mu = estimates(1);
eshock_std = estimates(2);
eshock_ro = estimates(3);
gain1 = estimates(4);
gain2 = estimates(5);

mu_tol = tol*mu;
eshock_std_tol = tol*eshock_std;
eshock_ro_tol = tol*eshock_ro;
gain1_tol = tol*gain1;
gain2_tol = tol*gain2;


%Calculate 'B'

B = zeros(5,5);

%Calculate numerical derivatives for mu
mu_ub_all = msm_g_st(seedstate,[mu+mu_tol,eshock_std,eshock_ro,gain1,gain2],sT,EM,expectationshock,constant_2,n);
mu_lb_all = msm_g_st(seedstate,[mu-mu_tol,eshock_std,eshock_ro,gain1,gain2],sT,EM,expectationshock,constant_2,n);
for count = 1:5
    B(count,1) = (mu_ub_all(count)-mu_lb_all(count))/(2*mu_tol);
end

%Calculate numerical derivatives for eshock_std
eshock_std_ub_all = msm_g_st(seedstate,[mu,eshock_std+eshock_std_tol,eshock_ro,gain1,gain2],sT,EM,expectationshock,constant_2,n);
eshock_std_lb_all = msm_g_st(seedstate,[mu,eshock_std-eshock_std_tol,eshock_ro,gain1,gain2],sT,EM,expectationshock,constant_2,n);
for count = 1:5
    B(count,2) = (eshock_std_ub_all(count)-eshock_std_lb_all(count))/(2*eshock_std_tol);
end

%Calculate numerical derivatives for eshock_ro
eshock_ro_ub_all = msm_g_st(seedstate,[mu,eshock_std,eshock_ro+eshock_ro_tol,gain1,gain2],sT,EM,expectationshock,constant_2,n);
eshock_ro_lb_all = msm_g_st(seedstate,[mu,eshock_std,eshock_ro-eshock_ro_tol,gain1,gain2],sT,EM,expectationshock,constant_2,n);
for count = 1:5
    B(count,3) = (eshock_ro_ub_all(count)-eshock_ro_lb_all(count))/(2*eshock_ro_tol);
end

%Calculate numerical derivatives for gain (agent 1)
gain1_ub_all = msm_g_st(seedstate,[mu,eshock_std,eshock_ro,gain1+gain1_tol,gain2],sT,EM,expectationshock,constant_2,n);
gain1_lb_all = msm_g_st(seedstate,[mu,eshock_std,eshock_ro,gain1-gain1_tol,gain2],sT,EM,expectationshock,constant_2,n);
for count = 1:5
    B(count,4) = (gain1_ub_all(count)-gain1_lb_all(count))/(2*gain1_tol);
end

%Calculate numerical derivatives for gain (agent 2)
gain2_ub_all = msm_g_st(seedstate,[mu,eshock_std,eshock_ro,gain1,gain2+gain2_tol],sT,EM,expectationshock,constant_2,n);
gain2_lb_all = msm_g_st(seedstate,[mu,eshock_std,eshock_ro,gain1,gain2-gain2_tol],sT,EM,expectationshock,constant_2,n);
for count = 1:5
    B(count,5) = (gain2_ub_all(count)-gain2_lb_all(count))/(2*gain2_tol);
end

end