% L_D_max.fcn calculates the maximum lift to drag ratio given the drag
% polar of the airplane. This is the maximum glide ratio for zero wind.
% Phillips, Mechanics of Flight, 2nd Ed., Equation 3.2.13.
%
% INPUTS:
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   L_D: Maximum lift to drag ratio
%
% Written by:
%   Sam Jaeger
%   2/1/2023

function L_D = L_D_max(C_D_0,C_D_1,C_D_2)
    L_D = sqrt(1/C_D_2)/(2*sqrt(C_D_0) + C_D_1*sqrt(1/C_D_2));
end