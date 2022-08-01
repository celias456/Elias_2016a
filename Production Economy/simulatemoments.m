function [stats_returns,stats_pdratio,stats_macro,stats_errors,stats_moments] = simulatemoments(input,T,totalsimulations,burn)
%Simulates the model and outputs simulated moments
%   Parameters (1x4) (input)
%   Number of time periods for each simulation T
%   Total number of simulations (totalsimulations)
%   Number of initial observations to remove in each simulation (burn)

% Misspecification
% Agent 2 omits the lag of equity price in its PLM

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


%DECLARATIONS
mu = input(1); %Proportion of agents in economy who have a correctly specified PLM (agent 1)
eshock_std = input(2); %Standard deviation of expectation shock
eshock_ro = input(3); %AR coefficient on expectation shock
gain = input(4); %Gain parameter for CG learning

ahoc = .9; %Value for initial condition of parameter estimates

eshock_1_on = 1; %Equals 1 if expectation shock for agent 1 is activated
eshock_2_on = 0; %Equals 1 if expectation shock for agent 2 is activated

seedstate = 7; %Set state value for pseudo random number generator
stream = RandStream('mt19937ar','Seed',seedstate);
RandStream.setDefaultStream(stream);

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

%Get Steady State Values
ss = get_steadystates(alpha,beta,delta,psi);

% Find REE solutions
[Pbar,Mbar1,Mbar2] = get_ree(a1,a2,b1,ro);

%Get gamma coefficients
[G] = get_gammas();
 
%% Initialize Variables for Simulation
%Starting values for parameter estimates; %First value is phi_c1, second value is phi_p, 
%third value is phi_z, fourth value is phi_c2, fifth value is theta
P = ahoc.*Pbar;  

data_returns = zeros(totalsimulations,6);
data_pdratio = zeros(totalsimulations,4);
data_macro = zeros(totalsimulations,10);
data_errors = zeros(totalsimulations,6);


for count_simulations = 1:totalsimulations
    
    %Initialize variables
    shock = sde*randn(T+2,1); %Shock values
    eshock_wn = eshock_std*randn(T+2,2); %Expectation shock white noise values
    P_0 = get_vmap(P,mu,a1,a2,b1,ro)*shock(1); %Initial value for equity price
    Z_0 = shock(1); %Initial value for shock 
    
    %Get law of motion for all variables (log deviation from the steady
    %state
    [allvariables,pf1_totaltimes,pf2_totaltimes,pf1_engaged,pf2_engaged,R1test,R2test] = lom_learn_cg(a1,a2,b1,ro,G,P,Mbar1,Mbar2,P_0,Z_0,shock,mu,T,gain,eshock,eshock_wn);
     
    %Get levels of all variables
    allvariables_levels = get_levels(ss,T,allvariables,burn);

    %Get statistics
    data_returns(totalsimulations,:) = statistics_returns(allvariables_levels);
    data_pdratio(totalsimulations,:) = statistics_pdratio(allvariables_levels);
    data_macro(totalsimulations,:) = statistics_macro(allvariables_levels);
    data_errors(totalsimulations,:) = [pf1_totaltimes,pf2_totaltimes,pf1_engaged,pf2_engaged,R1test,R2test];
end


%Calculate simulation statistics
[stats_returns,stats_pdratio,stats_macro,stats_errors,stats_moments] = statistics_simulation(totalsimulations,data_returns,data_pdratio,data_macro,data_errors);

end