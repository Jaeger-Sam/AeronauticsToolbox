% stall.fcn computes the stall speed given the maximum lift
% coefficient and the design ratio (f_1). The stall speed units will be
% given by the units of f_1 inputted by the user. For instance if the user
% uses lb for weight, ft^2 for area, and slugs/ft^3 for density the
% resulting stall speed will be in ft/s. Equation from Mechanics of Flight,
% Phillips, 2nd Ed., (3.8.3).
%
% INPUTS:
%   C_L_max: maximum vehicle lift coefficient (normalize by main wing area)
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density 
%
% OUTPUTS:
%   V_min: stall speed
%
% Written by:
%   Sam Jaeger
%   1/31/2023

function V_min = stall(C_L_max,f_1)
    V_min = (C_L_max^(-1/2))*sqrt(f_1); 
end