% Thin_Airfoil computes the zero angle of attack and moment coefficient at
% the quarter chord based off of a NACA 4 digit series airfoil.
%
% Inputs:
%   m: maximum camber in % of chord (camber_max = c*m/100)
%   p: location of maximum camber in tenths (x/c = p/10)
%   c: length of chord
%
% Outputs:
%   alpha_L0: zero lift angle of attack (degrees)
%   C_m_c4: section moment coefficent at the quarter chord
%
% Written by:
%   Sam Jaeger
%   Revised: 11/11/2022

function [alpha_L0, C_m_c4] = Thin_Airfoil(m,p,c)

    % point to break integral up in theta coords
    theta_1 = acos( 1 - 2*(p/10)/c);
  
    %% alpha_L0
    int_1 = integral(@(theta) NACA_4_dig_dyc_dx_1(theta,m,p),0,theta_1);
    int_2 = integral(@(theta) NACA_4_dig_dyc_dx_2(theta,m,p),theta_1,pi);
    
    alpha_L0 = (int_1 + int_2)/pi;
    alpha_L0 = alpha_L0*180/pi;
    
    %% C_m_c_4
    
    int_3 = integral(@(theta) NACA_4_dig_dyc_dx_3(theta,m,p),0,theta_1);
    int_4 = integral(@(theta) NACA_4_dig_dyc_dx_4(theta,m,p),theta_1,pi);
    
    C_m_c4 = (int_3 + int_4)*0.5;
end