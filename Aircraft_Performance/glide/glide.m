% glide.fcn calculates the best glide airspeed and associated glide ratio
% with wind condition and airplane performance. This assumes small glide 
% angles. User must keep track of units. 
% From Phillips, Mechanics of Flight, 2nd ed., Equation 3.7.22 and 3.7.21.
%
% INPUTS:
%   V_w: wind speed
%   phi_w: wind track angle (rad). Angle between the ground track and the
%           wind direction.
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho = air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   V_g: best glide airspeed
%   R_g: glide ratio (== ground speed / sink rate)
%
% Written by:
%   Sam Jaeger
%   2/1/2023

function [V_g, R_g] = glide(V_w,phi_w,f_1,C_D_0,C_D_1,C_D_2)
    a = C_D_1/C_D_2;
    b = 12*C_D_0/C_D_2;

    V_MDV = 2/sqrt(2)*f_1/sqrt(a + sqrt(a^2 + b));

    fun = @(V)(2*((V^4)-(V_MDV^4))*V - ((3*V^4)-(V_MDV^4))*V_w*(cos(phi_w)*sqrt(1- ((V_w*sin(phi_w)^2)/(V^2) )) + (V_w/V*sin(phi_w)^2) ) );
    V_g = fzero(fun, V_MDV);

    num = sqrt(1 - (((V_w*sin(phi_w))^2)/(V_g^2)) ) - V_w*cos(phi_w)/V_g;
    denom = (V_g^2)*C_D_0/f_1 + C_D_1 + f_1*C_D_2/(V_g^2);

    R_g = num / denom;
end