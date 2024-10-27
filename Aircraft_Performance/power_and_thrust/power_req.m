% Power_req.fcn calculates the power required for a certain flight
% condition and aircraft. The user must keep track of units. From Mechanics
% of Flight, Phillips, 2nd Ed., Equation (3.3.7).
% 
% INPUTS:
%   V: velocity
%   W: weight
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   P_R: Power required
%
% Written by:
%   Sam Jaeger
%   1/31/2023

function P_R = power_req(V,W,f_1,C_D_0,C_D_1,C_D_2)
   P_R = (C_D_0/f_1*(V^3) + C_D_1*V + f_1*C_D_2/V )*W;
end