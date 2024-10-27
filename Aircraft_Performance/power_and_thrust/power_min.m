% Power_min.fcn calculates the minimum power for a set of aircraft 
% parameters. The user must keep track of units. From Mechanics
% of Flight, Phillips, 2nd Ed., Equation (3.3.6), (3.3.10).
% 
% INPUTS:
%   W: weight
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   P_R_min: Minimum Power required
%
% Written by:
%   Sam Jaeger
%   1/31/2023

function P_R_min = power_min(W,f_1,C_D_0,C_D_1,C_D_2)
   C_L = (1/2/C_D_2)*(C_D_1 + sqrt(C_D_1^2 + 12*C_D_0*C_D_2));
   P_R_min = W*sqrt(f_1)*( (C_D_0/(C_L^(3/2))) + (C_D_1/(C_L^(1/2))) + (C_D_2*(C_L^(1/2))));
end