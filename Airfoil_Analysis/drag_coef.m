% drag_coef.fcn computes the drag coefficient given a lift coefficient and
% the corresponding drag coefficient constants.
%   C_D = C_D_0 + C_D_1*C_L + C_D_2*C_L^2
%
% INPUTS:
%   C_L: vector of angles of attacks (rad)
%   C_D_0: constant drag coef.
%   C_D_1: linear drag coef.
%   C_D_2: parobolic drag coef.
%
% INPUTS
%   C_D: vector of drag coefficients corresponding to angles of attack
%
% Sam Jaeger
% 1/15/2024

function C_D = drag_coef(C_L, C_D_0, C_D_1, C_D_2)
    C_D = C_D_0 + C_D_1.*C_L + C_D_2.*(C_L.^2);
end