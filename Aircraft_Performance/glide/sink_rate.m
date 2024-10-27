% sink_rate.fcn calculates the sink rate for a certain flight
% condition and aircraft. The user must keep track of units. From Mechanics
% of Flight, Phillips, 2nd Ed., Equation (3.7.3).
% 
% INPUTS:
%   V: velocity
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   V_s: sink rate
%
% Written by:
%   Sam Jaeger
%   2/1/2023

function V_s = sink_rate(V,f_1,C_D_0,C_D_1,C_D_2)
    V_s = C_D_0*(V^3)/f_1 + C_D_1*V + f_1*C_D_2;
end