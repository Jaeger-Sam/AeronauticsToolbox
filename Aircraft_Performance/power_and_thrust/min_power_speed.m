% min_power_speed.fcn calculates the minimum power airspeed given 
% parameters of the aircraft. This is the same speed that minimizes sink
% rate.
% The user must keep track of units. 
% From Mechanics of Flight, Phillips, 2nd Ed., Equation (3.3.12).
% 
% INPUTS:
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   V_MDV: Minimum power velocity
%
% Written by:
%   Sam Jaeger
%   1/31/2023

function [V_MDV] = min_power_speed(f_1,C_D_0,C_D_1,C_D_2)
    a = C_D_1/C_D_2;
    b = 12*C_D_0/C_D_2;

    V_MDV = 2/sqrt(2)*f_1/sqrt(a + sqrt(a^2 + b));
end