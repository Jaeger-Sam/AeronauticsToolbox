% climb_rate.fcn calclates the climb rate given flight conditions and 
% performance parameters. The user must keep track of units. From Mechanics
% of Flight, Phillips, 2nd Ed., section 3.4. Must assume small climb angle!
% 
% INPUTS:
%   V: airspeed
%   PW: power loading (equiv P/W)
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
%   Revised: 2/6/23
%       Small angle assumption
%   Revised 2/13/23
%       Power loading 

function V_c = climb_rate(V,PW,f_1,C_D_0,C_D_1,C_D_2)
    V_c = PW - (C_D_0*(V^3)/f_1 + C_D_1*V + C_D_2*f_1/V);
end