function [g_sT] = g_st(input,sT,EM)
%Simulates the model and outputs deviations of the simulated moments from
%the empirical moments

%   Parameters (1x4) input
%   Number of Simulations sT
%   Empirical moments (1x4) EM

% Misspecification
% Agent 2 omits the productivity shock in its PLM

% Expectation shock is AR(1) process 

% Variables   1. equity price 
%             2. output         
%             3. consumption    
%             4. dividends      
%             5. bond price     
%             6. capital
%             8. investment
%             9. gross profit
%
% State Variables    shock          
%                    equity price  

% Funtions used by this code
% STEADYSTATES.m
% SOLVE.m
% i.  TMAP.m
% ii. VMAP.m
% GAMMAS.m
% VMAP.m
% LOM_HEE.m
% LEARN_LS.m
% LEARN_CG.m
% LEVELS.m
% STATS.m
% i. RETURNS.m
% SIMSTATS.m


%DECLARATIONS
mu = input(1); %Proportion of agents in economy who have a correctly specified PLM (agent 1)
eshock_std = input(2); %Standard deviation of expectation shock
eshock_ro = input(3); %AR coefficient on expectation shock
gain = input(4); %Gain parameter for CG learning

totalsimulations = 1; %Number of simulations to run
eshock_1_on = 1; %Equals 1 if expectation shock for agent 1 is activated
eshock_2_on = 0; %Equals 1 if expectation shock for agent 2 is activated
burn = sT/2; %Number of initial observations to remove

seed_state = 98; %Set state value for pseudo random number generator
randn('state',seed_state);
ahoc = .9; %Value for Ad-hoc initial condition of Phi AHOC-B

%Deep parameters of the model (from Mathematica notebook)
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

%Expectation Shock parameters
eshock_1_mean = 0; %Mean of expectation shock for agent 1
eshock_1_std = eshock_std; %Standard deviation of expectation shock for agent 1
eshock_2_mean = 0; %Mean of expectation shock for agent 2 
eshock_2_std = eshock_std;  %Standard deviation of expectation shock for agent 2
eshock = [eshock_1_on,eshock_1_mean,eshock_1_std,eshock_ro;eshock_2_on,eshock_2_mean,eshock_2_std,eshock_ro];

%GET STEADY STATE VALUES
ss = STEADYSTATES(alpha,beta,delta,psi);

%FIND REE/HEE SOLUTIONS
[Pbar,Mbar1,Mbar2] = SOLVE(a1,a2,b1,ro);

%GET GAMMA COEFFICIENTS
[G] = GAMMAS();
 
%INITIALIZE VARIABLES FOR SIMULATION
%Starting values for parameter estimates; %First value is phi_c1, second value is phi_p, 
%third value is phi_z, fourth value is phi_c2, fifth value is theta
P = ahoc*Pbar; 

data = zeros(totalsimulations,33);
moments_simulated = zeros(totalsimulations,5);


for count_simulations = 1:totalsimulations
    
    %Initialize variables
    shock = sde*randn(sT+2,1); %Shock values
    eshock_wn = eshock_std*randn(sT+2,2); %Expectation shock white noise values
    P_0 = VMAP(P,mu,a1,a2,b1,ro)*shock(1); %Initial value for equity price
    Z_0 = shock(1); %Initial value for shock 

    %GET LAW OF MOTION FOR ALL VARIABLES (LOG DEVIATIONS FROM THE SS)
    [allvariables,pf1_totaltimes,pf2_totaltimes,pf1_engaged,pf2_engaged,R1test,R2test] = LEARN_CG(a1,a2,b1,ro,G,P,Mbar1,Mbar2,P_0,Z_0,shock,mu,sT,gain,eshock,eshock_wn);
     
    %GET LEVELS OF VARIABLES
    allvariables_levels = LEVELS(ss,sT,allvariables,burn);

    %GET STATISTICS
    [s,simmoments] = STATS(allvariables_levels);
    data(count_simulations,:) = [s,pf1_totaltimes,pf2_totaltimes,pf1_engaged,pf2_engaged,R1test,R2test];
    moments_simulated(count_simulations,:) = simmoments;
end


%CALCULATE SIMULATION STATISTICS
stats = SIMSTATS(data,moments_simulated);

%Calculate criterion function
g_sT = [stats.mom_1-EM(1);stats.mom_2-EM(2);stats.mom_3-EM(3);stats.mom_4-EM(4)];


end

