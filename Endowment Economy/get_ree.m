function [Pbar,Mbar1,Mbar2] = get_ree(ro,beta,gamma,constant_2) 
% Returns HEE solutions
% Pbar is a 2x1 vector
%     P(1) = REE constant for agent 1
%     P(2) = REE coefficient on dividends
%
% Mbar1 is a 2x2 matrix used as an initial value for R1 in AL
% Mbar2 is either a scalar or 2x2 matrix used as the initial value for the
% R2 matrix in AL
%
% -------------------------------------------------------------------------

% Find solution
phibar_c = 0;
phibar_d = ((ro*(1-beta-gamma)+gamma)/(1-beta*ro)); %coefficient on d
        
% Pbar are the solutions to be learned with phi_c1 listed first and phi_d
% listed second

Pbar = [phibar_c;phibar_d];

%Mbar1 and Mbar2 are used as the initial values for the R matrices in AL
Mbar1 = .01.*ones(2,2);

if constant_2 == 1
    Mbar2 = .01.*ones(2,2);
else
    Mbar2 = .01;
end

end 
