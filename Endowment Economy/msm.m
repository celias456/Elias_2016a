% Endowment based asset pricing model

% Simulated Method of Moments Estimation

% Learning Misspecification
%   Agent 2 includes lagged equity price in its PLM
%   Agent 2 does not include a constant in it's PLM

%Order of moments
%1. Equity Return Mean
%2. Equity Return Variance
%3. Price-Dividend Ratio Mean
%4. Price-Dividend Ratio Variance
%5. Price-Dividend Ratio Autocovariance

%Order of parameters
%1. Degree of heterogeneity
%2. Standard deviation of expectation shock
%3. AR coefficient on expectation shock
%4. learning gain parameter (agent 1)
%5. learning gain parameter (agent 2)

% =========================================================================
clear all
clc;
warning('off','all');

%% Set specifications for estimation
T = 256; %Number of time periods in empirical data
s = 5; %Multiple of T periods to simulate (should be greater than one)
n = 2000; %Number of replications 
tol = .00001; %tolerance percentage for calculation of standard errors

%Initial values used in estimation
mu_initial = 1; 
eshock_std_initial = .015;
esock_ro_initial = .95;
gain1_initial = .3;
gain2_initial = 0;

%Indicator that turns on expectation shock
%1 = on, 0 = off
%first element is for agent 1, second element is for agent 2
expectationshock = [1,0]; 

%Indicator of type of weighting matrix to use
%1 = identity matrix
%2 = optimal matrix
weightmatrix_type = 1;

%Indicator for a constant in agent type two's PLM
%1 = constant
%0 = no constant
constant_2 = 1;

%Set seed state value for pseudo random number generator
seedstate = 11; 

%Calculation of sT (should be integer greater than one)
sT = s*T; %Number of time periods in each simulation

%Get empirical moments observed in data
em1 = 2.0145; %equity return mean
em2 = 65.7203; %equity return variance
em3 = 34.1853; %P-D ratio mean
em4 = 257.1409; %P-D ratio variance
em5 = 251.6735; %P-D ratio autocovariance
EM = [em1,em2,em3,em4,em5];

%Calculate weighting matrix W
[W,omega_hat] = msm_weightmatrix(weightmatrix_type);

%% Specify bounds and initial values for criterion function minimization
mu_lb = .8;
mu_ub = .95;

eshock1_std_lb = .01;
eshock1_std_ub = .02;

eshock2_std_lb = .01;
eshock2_std_ub = .02;

eshock1_ro_lb = .94;
eshock1_ro_ub = .99;

eshock2_ro_lb = .92;
eshock2_ro_ub = .97;

gain1_lb = .01;
gain1_ub = .5;

gain2_lb = .01;
gain2_ub = .5;

theta_initial = [eshock_std_initial,esock_ro_initial,gain1_initial];
[cfvalue_theta_initial,g_st_theta_initial,test] = msm_cf(seedstate,theta_initial,sT,EM,W,expectationshock,constant_2,n);

%% Estimate parameters
%Constraints/Options on minimization
A = [];
b = [];
Aeq = [];
beq = [];
if expectationshock(1) == 1
    lb = [eshock1_std_lb,eshock1_ro_lb,gain1_lb];
    ub = [eshock1_std_ub,eshock1_ro_ub,gain1_ub];
end
if expectationshock(2) == 1
    lb = [mu_lb,eshock2_std_lb,eshock2_ro_lb,gain1_lb,gain2_lb];
    ub = [mu_ub,eshock2_std_ub,eshock2_ro_ub,gain1_ub,gain2_lb];
end 
options = optimset('Algorithm','interior-point');
%Assign function handle
f = @(input)msm_cf(seedstate,input,sT,EM,W,expectationshock,constant_2,n);
%Perform minimization
[theta_hat,cfvalue_theta_hat] = fmincon(f,theta_initial,A,b,Aeq,beq,lb,ub,[],options);


%% Calculate Standard errors

% B = msm_numderiv(seedstate,theta_hat,tol,sT,EM,expectationshock,constant_2,n);
% cov1 = ((B'*W*B)^-1);
% theta_hat_cov = cov1*B'*W*(1+(1/s))*omega_hat*W*B*cov1';
% theta_hat_se = [sqrt(theta_hat_cov(1,1)),sqrt(theta_hat_cov(2,2)),sqrt(theta_hat_cov(3,3)),sqrt(theta_hat_cov(4,4)),sqrt(theta_hat_cov(5,5))];


%% Simulate model

% [stats_returns,stats_pdratio,stats_predictability,stats_macro,stats_errors,stats_moments,stats_prediction_errors] = msm_model(n,sT,sT/2,theta_hat,seedstate);

