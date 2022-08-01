function [gamma] = get_gammas() 
%Obtains gamma coefficients for endogeneous variables in model
%gamma is a 5x4 matrix holding the gamma coefficients for the model
%    gamma(1,:) = gamma coefficients for output
%    gamma(2,:) = gamma coefficients for consumption
%    gamma(3,:) = gamma coefficients for dividends
%    gamma(4,:) = gamma coefficients for bond price
%    gamma(5,:) = gamma coefficients for capital
%    gamma(6,:) = gamma coefficients for investment
%    gamma(7,:) = gamma coefficients for gross profit


load gammas_pe.txt

gamma = gammas_pe; %coefficients for law of motion of endogeneous variables in the model


end