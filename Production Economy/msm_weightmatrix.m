function [W,omega_hat] = msm_weightmatrix(weightmatrix_type)
%Calculates the weight matrix W

load NW.csv
NW_data = NW;
T = size(NW_data,1);
p = round((T^.25));
omega_hat = covnw(NW_data,p,0);

if weightmatrix_type == 1
    W = eye(3);
else
    W = omega_hat^-1;

end

end