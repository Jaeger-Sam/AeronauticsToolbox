% min_drag_speed.fcn calculates the minimum drag airspeed for a certain flight
% condition and aircraft. This is the velocity where L/D is maximized and
% is the best glide speed with zero wind.
% The user must keep track of units. From Mechanics
% of Flight, Phillips, 2nd Ed., Equation (3.2.14).
% 
% INPUTS:
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%   alpha_T: thrust mounting angle (with respect to freestream velocity)
%
% OUTPUS:
%   V_MD: minimum-drag velocity
%
% Written by:
%   Sam Jaeger
%   2/1/2023

function V_MD = min_drag_speed(f_1,C_D_0,C_D_1,C_D_2,alpha_T)
    V_MD = ((C_D_0/C_D_2)^(-1/4))*sqrt(f_1)*sqrt(1/(1 + ((2*sqrt(C_D_0*C_D_2)+C_D_1)*tan(alpha_T)) ));
end