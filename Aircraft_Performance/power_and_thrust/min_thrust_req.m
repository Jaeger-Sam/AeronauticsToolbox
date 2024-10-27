% min_thrust_req.fcn calculates the minimum thrust required for a certain flight
% condition and aircraft. This is the velocity where L/D is maximized and
% at V_MD. Thrust angle alpha_T must be small.
% The user must keep track of units. From Mechanics
% of Flight, Phillips, 2nd Ed., Equation (3.2.26).
% 
% INPUTS:
%   W: weight of the aircraft
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   T_R_min: minimum thrust required
%
% Written by:
%   Sam Jaeger
%   2/1/2023

function T_R_min = min_thrust_req(W,C_D_0,C_D_1,C_D_2)
    T_R_min = W*(2*sqrt(C_D_0*C_D_2) + C_D_1);
end