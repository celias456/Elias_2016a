function [Pbar,Mbar1,Mbar2] = get_ree(a1,a2,b1,ro,constant_2) 
% Returns REE solutions
% Pbar is a 5x1 matrix
%     P(1) = REE constant (agent 1)
%     P(2) = REE coefficient on p (equity price) for agent 1
%     P(3) = REE coefficient on z (productivity shock) for agent 1
%     P(4) = REE constant (agent 2)
%     P(5) = REE coefficient on p (equity price) for agent 2
%
% Mbar1 is a 3x3 matrix used as the initial values for the R1 matrix in AL
% Mbar2 is either a scalar or 2x2 matrix used as the initial value for the R2 matrix in AL
%  
% -------------------------------------------------------------------------


% Find REE solutions
phibar_p = (1 - sqrt(1-4*a1*a2))/(2*a1); %coefficient on price
phibar_z = (b1*ro)/(1-a1*(ro+((1 - sqrt(1-4*a1*a2))/(2*a1)))); %coefficient on productivity shock

thetabar = phibar_p; %This is for agents of type 2

% Pbar are the solutions to be learned with phi_c1 listed first, phi_p
% listed second, phi_z listed third, phi_c2 listed fourth, and theta listed
% fifth
Pbar = [0;phibar_p;phibar_z;0;thetabar];

%Mbar1 and Mbar2 are used as the initial values for the R matrices in AL
Mbar1 = .01.*ones(3,3);
if constant_2 == 1
    Mbar2 = .01.*ones(2,2);
else
    Mbar2 = .01;
end

end   