% speed.fcn calculates the speed given power and condition of the 
% aircraft. The minimum power velocity is used as a first guess to solve 
% for the final velocity. The user must keep track of units. 
% From Mechanics of Flight, Phillips, 2nd Ed., Equation (3.3.7), (3.3.12).
% 
% INPUTS:
%   PW: Power loading (equiv P/W)
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   V_inf: Speed
%   V_MDV: Minimum power velocity
%
% Written by:
%   Sam Jaeger
%   1/31/2023
%   Revised: 2/13/2023
%       Power loading

function [V_inf,V_MDV] = speed(PW,f_1,C_D_0,C_D_1,C_D_2)
    a = C_D_1/C_D_2;
    b = 12*C_D_0/C_D_2;

    V_MDV = 2/sqrt(2)*f_1/sqrt(a + sqrt(a^2 + b));

    V_r = roots([C_D_0/f_1/PW, 0, C_D_1/PW, -1, f_1*C_D_2/PW]);
    V_inf = V_r(3);
%   V_inf = fzero(@(V)( (C_D_0/f_1*(V^3) + C_D_1*V + f_1*C_D_2/V )*((PW)^(-1)) - 1),V_MDV*10);
   
%    if V_inf < 0
%        PW
%        error('Negative Velocity! Check units!')
%    end
end