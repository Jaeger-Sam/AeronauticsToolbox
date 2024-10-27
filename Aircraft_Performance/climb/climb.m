% climb.fcn calclates the maximum rate of climb, rate of climb at the
% optimal angle of climb, as well as Vy and Vx. Vy and Vx are first
% calculated via a root finder. Vy is found by taking a derivative and then
% setting it to zero. Vx is found by finding the point where a tangent line
% orginating from the orgin intersects with the rate of climb curve. Then 
% the associated climb rates are found via climb_rate.fcn. 
% climb_rate.fcn must be in the same folder as this function.
% The user must keep track of units! From Mechanics of Flight, Phillips, 
% 2nd Ed., section 3.4. Additionally, see: 
% https://www.kitplanes.com/using-level-accelerations-to-determine-climb-performance/
% 
% Must assume small climb angle!
% 
% INPUTS:
%   PW: power loading (equiv P/W)
%   f_1: design variable. f_1 == 2*(W/S_w)/rho where W = vehicle weight,
%           S_w == main wing area, rho == air density
%   C_D_0: Parasidic drag coef
%   C_D_1: Linear drag coef
%   C_D_2: Induced drag coef. C_D_2 == 1/(pi*e*R_A)
%
% OUTPUS:
%   V_c_max: maximum rate of climb
%   V_y: speed for maximum rate of climb
%   V_c_ang: rate of climb for best angle of climb
%   V_x: Best angle of climb
%
% Written by:
%   Sam Jaeger
%   2/6/23
%   Revised 2/13/23
%       Power loading 


function [V_c_max, V_y, V_c_ang, V_x] = climb(PW,f_1,C_D_0,C_D_1,C_D_2)
    Vy_roots = roots([-3*C_D_0/f_1, 0, -C_D_1, 0, C_D_2*f_1]);
    Vx_roots = roots([2*C_D_0/f_1, 0, 0, PW, -2*C_D_2*f_1]);
    
    V_y = Vy_roots(4);

    V_c_max = climb_rate(V_y,PW,f_1,C_D_0,C_D_1,C_D_2);

    V_x = Vx_roots(4);

    V_c_ang = climb_rate(V_x,PW,f_1,C_D_0,C_D_1,C_D_2);
end