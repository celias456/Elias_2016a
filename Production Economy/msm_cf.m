function [CF,g_sT,test] = msm_cf(seedstate,input,sT,EM,W,expectationshock,constant_2,n)
%Simulates the model and outputs criterion function in msm estimation
%   number for seedstate (seedstate)
%   Parameters (1x3) (input)
%   Number of time periods for each simulation (sT)
%   Empirical moments (1x5) (EM)
%   Weighting Matrix (W) 
%   Indicator for turning on expectation shocks (1x2) (expectationshock)
%   Indicator for constant in agent type 2's PLM (constant_2)
%   Number of replications (n)

stream = RandStream('mt19937ar','Seed',seedstate);
% RandStream.setDefaultStream(stream);
RandStream.setGlobalStream(stream);

%% Declarations
mu = 1; %Proportion of agents in economy who have a correctly specified PLM (agent 1)
eshock_std = input(1); %Standard deviation of expectation shock
eshock_ro = input(2); %AR coefficient on expectation shock
gain1 = input(3); %Gain parameter for CG learning (agent 1)
gain2 = 0; %Gain parameter for CG learning (agent 2)
gain = [gain1,gain2];

eshock_1_on = expectationshock(1,1); %Equals 1 if expectation shock for agent 1 is activated
eshock_2_on = expectationshock(1,2); %Equals 1 if expectation shock for agent 2 is activated

burn = sT-256; %Number of initial observations to remove

%% Deep parameters of the model (from Mathematica notebook)
ro  = 0.95; %AR parameter on shock
sde = 0.00712; %Standard deviation of the shock term
alpha = 0.36;
beta = 0.99;
delta = 0.025;
psi = (1-beta+delta*beta)/(alpha*beta);

%Coefficients of the reduced form  k(t) = a1*E(k(t+1)+a2*k(t-1)+b1*z(t)
a1 = 0.497089;    
a2 = 0.50211;    
b1 = 0.00361316;

%% Expectation Shock parameters
eshock_1_mean = 0; %Mean of expectation shock for agent 1
eshock_1_std = eshock_std; %Standard deviation of expectation shock for agent 1
eshock_2_mean = 0; %Mean of expectation shock for agent 2 
eshock_2_std = eshock_std;  %Standard deviation of expectation shock for agent 2
eshock = [eshock_1_on,eshock_1_mean,eshock_1_std,eshock_ro;eshock_2_on,eshock_2_mean,eshock_2_std,eshock_ro];

%% Get steady state values
ss = get_steadystates(alpha,beta,delta,psi);

%% Get ree solutions
[Pbar,Mbar1,Mbar2] = get_ree(a1,a2,b1,ro,constant_2);

%% Get gamma coefficients
[G] = get_gammas();

%% Initialize Variables for Simulation
 
%Starting values for parameter estimates; 
%First value is phi_c1, second value is phi_p, third value is phi_z, fourth value is phi_c2, fifth value is theta
ahoc = .9; %Value for initial condition of parameter estimates
P = ahoc*Pbar; 

%Create storage for simulation runs
r_mean = zeros(n,1);
r_var = zeros(n,1);
pd_mean = zeros(n,1);
pd_var = zeros(n,1);
pd_auto = zeros(n,1);
dc_var = zeros(n,1);
dd_var = zeros(n,1);

%% Begin Simulation

for i = 1:n
    
    shock = sde*randn(sT+2,1); %Productivity shock values
    eshock_wn = eshock_std*randn(sT+2,2); %Expectation shock white noise values
    P_0 = get_vmap(P,mu,a1,a2,b1,ro)*shock(1); %Initial value for equity price
    Z_0 = shock(1); %Initial value for shock

    %Get law of motion for all variables (log deviation from the ss)
    [allvariables] = lom_learn_cg(a1,a2,b1,ro,G,P,Mbar1,Mbar2,P_0,Z_0,shock,mu,sT,gain,eshock,eshock_wn,constant_2);

    %Get levels of all variables
    allvariables_levels = get_levels(ss,sT,allvariables,burn);

    %Get Moments
    [r_mean(i),r_var(i),pd_mean(i),pd_var(i),pd_auto(i),dc_var(i),dd_var(i)] = msm_calculate_moments(allvariables_levels);

end


%% Calculate criterion function

r1 = median(r_mean);
r2 = median(r_var);
pd1 = median(pd_mean);
pd2 = median(pd_var);
pd3 = median(pd_auto);
dc = median(dc_var);
dd = median(dd_var);

test = [r_mean,r_var,pd_mean,pd_var,pd_auto,dc_var,dd_var];

mom_diff_1 = r1-EM(1);
mom_diff_2 = r2-EM(2);
mom_diff_3 = pd1-EM(3);
mom_diff_4 = pd2-EM(4);
mom_diff_5 = pd3-EM(5);

g_sT = [mom_diff_1;mom_diff_3;mom_diff_4];

CF = (g_sT'*W*g_sT);

end