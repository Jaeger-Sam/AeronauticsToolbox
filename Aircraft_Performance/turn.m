% turn.fcn calculates the maneuvering speed, maximum turning rate, and
% maximum bank angle obtainable. Units must be in SI and be kept track by
% the user. From Mechanics of Flight, Phillips, 2nd Edition, Equations
% 3.9.17, 3.9.25, 3.9.33.
%
% INPUTS:
%   C_L_max: maximum lift coefficient in turn
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   f_1_max: design variable. f_1 == 2*(W_max/S_w)/rho where W_max = 
%               maximum vehicle weight, S_w = main wing area, rho=air density
%   n_pll: positive load factor
%   W: weight
%   W_max: maximum possible weight
%
% OUTPUTS:
%   V_M: maneuvering speed (maximum turns can be performed (also known as
%       the corner velocity and is the same velocity that gives the maximum
%       turning rate
%   Omega_max: maximum turning rate (rad/s)
%   phi_max: load limited bank angle (rad)
%
% Written by:
%   Sam Jaeger
%   2/1/2023

function [V_M,Omega_max,phi_max] = turn(C_L_max, f_1, f_1_max, n_pll, W, W_max)
    g = 9.802; %m/s^2

    V_M = sqrt(n_pll/C_L_max)*sqrt(f_1_max);
    Omega_max = sqrt(C_L_max*( (n_pll*W_max/W) - (W/n_pll/W_max) ) )*g*sqrt(1/f_1);
    phi_max = acos(W/n_pll/W_max);

end