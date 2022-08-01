% Endowment based asset pricing model

% Model Simulation

% Learning Misspecification
%   Agent 2 includes lagged equity price in its PLM

% Expectation shock is AR(1) process 

% Variable order
%1. equity price 
%2. output         
%3. consumption    
%4. dividends      
%5. bond price     
%6. capital  
%7. productivity shock
%8. investment
%9. gross profit

%State Variables    
%shock          
%equity price  


clear all
clc;
warning('off','all');


%% Simulation Specifications

n = 1; % Number of simulations
T = 11003;
burn = 1000;
% T = 256*5; %Number of time periods in each simulation
% burn = T-256; %Number of initial observations to remove


%LOM indicator: 1 = REE; 2 = LS Learning; 3 = CG Learning
lom_type = 1;

eshock_1_on = 0; %Equals 1 if expectation shock for agent 1 is activated
eshock_2_on = 0; %Equals 1 if expectation shock for agent 2 is activated

constant_2 = 1; %Equals 1 if agent type 2 has a constant in PLM, 0 otherwise

%% Parameters of Interest
%mu = 1;
mu = 1; %Proportion of agents in economy who have a correctly specified PLM (agent 1)

eshock1_std = .0166; %Standard deviation of the expectation shock for agent 1
eshock1_ro = .95; %AR coefficient on expectation shock for agent 1

eshock2_std = 0; %Standard deviation of the expectation shock for agent 2 
eshock2_ro = 0; %AR coefficient on expectation shock for agent 2 

gain1 = .3196; %Gain parameter for CG learning for agent 1
gain2 = .0100; %Gain parameter for CG learning for agent 2
gain = [gain1,gain2];

%% Fix random number generator
seedstate = 11; %Set state value for pseudo random number generator
stream = RandStream('mt19937ar','Seed',seedstate);
%RandStream.setDefaultStream(stream);
RandStream.setGlobalStream(stream);

%% Deep parameters of the model (from Mathematica notebook)
ro  = 0.95; %AR parameter on shock
sde = 0.00712; %Standard deviation of the shock term
alpha = 0.36;
beta = 0.99;
delta = 0;
gamma = 1;
psi = (1-beta+delta*beta)/(alpha*beta);
%Coefficients of the reduced form  p(t) = a*E(p(t+1)+b*d(t)
a = beta;        
b = (1-beta-gamma)*ro + gamma;

%% Expectation Shock parameters
eshock_1_mean = 0; %Mean of expectation shock for agent 1
eshock_1_std = eshock1_std; %Standard deviation of expectation shock for agent 1
eshock_2_mean = 0; %Mean of expectation shock for agent 2 
eshock_2_std = eshock2_std;  %Standard deviation of expectation shock for agent 2
eshock_1_ro = eshock1_ro; %First-order autocorrelation term for expectation shock for agent 1
eshock_2_ro = eshock2_ro; %First-order autocorrelation term for expectation shock for agent 2
eshock = [eshock_1_on,eshock_1_mean,eshock_1_std,eshock_1_ro;eshock_2_on,eshock_2_mean,eshock_2_std,eshock_2_ro];

%% Get Steady State Values
ss = get_steadystates(alpha,beta,delta,psi);

%% Find REE solutions
[Pbar,Mbar1,Mbar2] = get_ree(ro,beta,gamma,constant_2);

%% Get Gamma Coefficients
[G] = get_gammas();
 
%% Initialize Variables for Simulation
ahoc = .01; %Value for initial condition of parameter estimates
P = ahoc.*Pbar; %Starting values for parameter estimates; %First value is phi_c1, second value is phi_d

%Create storage for simulation runs
data_returns = zeros(n,7);
data_pdratio = zeros(n,5);
data_predictability = zeros(n,9);
data_macro = zeros(n,10);
data_errors = zeros(n,4);
data_predict_errors = zeros(n,2);

for i = 1:n
    
    %Generate productivity shocks and expectation shocks
    shock = sde*randn(T+2,1); %Shock values
    eshock_1_wn = eshock_1_std*randn(T+2,1); %Expectation shock white noise values for agent 1
    eshock_2_wn = eshock_2_std*randn(T+2,1); %Expectation shock white noise values for agent 2
    eshock_wn = [eshock_1_wn,eshock_2_wn]; %Expectation shock white noise values
    
    %Calculate initial values for equity price and productivity shock
    P_0 = get_vmap(P,mu,a,b,ro)*shock(1); %Initial value for equity price
    Z_0 = shock(1); %Initial value for shock

    %Get law of motion for all variables (log deviation from the steady
    %state
    if lom_type == 1
        [allvariables,pf_totaltimes,pf_engaged,R1test,R2test] = lom_ree(a,b,ro,G,Pbar,Z_0,shock,mu,T,eshock,eshock_wn);
    elseif lom_type == 2
        [allvariables,pf_totaltimes,pf_engaged,R1test,R2test] = lom_learn_ls(a,b,ro,G,P,Mbar1,Mbar2,P_0,Z_0,shock,mu,T,eshock,eshock_wn,constant_2);
    else
        [allvariables,pf_totaltimes,pf_engaged,R1test,R2test,Phi,Theta] = lom_learn_cg(a,b,ro,G,P,Mbar1,Mbar2,P_0,Z_0,shock,mu,T,gain,eshock,eshock_wn,constant_2);
    end
    
    %Get levels of variables
    allvariables_levels = get_levels(ss,T,allvariables,burn);

    %Get Statistics
    %[data_returns(i,:),r_gross,rf_gross] = statistics_returns(allvariables_levels);
    [data_pdratio(i,:),pdratio] = statistics_pdratio(allvariables_levels);
    %data_predictability(i,:) = statistics_predictability(r_gross,rf_gross,pdratio);
    %data_macro(i,:) = statistics_macro(allvariables_levels);
    %data_errors(i,:) = [pf_totaltimes,pf_engaged,R1test,R2test];
    %data_predict_errors(i,:) = statistics_predictions(allvariables_levels);
    
end


%% Calculate Simulation Statistics
%[stats_returns,stats_pdratio,stats_predictability,stats_macro,stats_errors,stats_moments,stats_prediction_errors] = statistics_simulation(data_returns,data_pdratio,data_predictability,data_macro,data_errors,data_predict_errors);


%% Learning Coefficients Graph

% Period = (1:10000);
% Constant = Phi(1,1002:end);
% ProdShock = Phi(2,1002:end);
% 
% 
% figure
% subplot(2,1,1)
% plot(Period,Constant)
% xlabel('Period')
% ylabel('Constant')
% 
% subplot(2,1,2)
% plot(Period,ProdShock)
% xlabel('Period')
% ylabel('Productivity Shock')