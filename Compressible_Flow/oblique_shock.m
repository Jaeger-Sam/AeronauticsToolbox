% [M2,Mn2,P0_rat,T_rat,P_rat,rho_rat] = oblique_shock(theta,M1,gamma)
%
% oblique_shock.m calculates the oblique shock properties of a perfect gas
% given a turning angle theta, upstream mach, and gamma. This function
% assumes the weak shock solution and calls the functions: theta_beta_M,
% normal_shock.m
%
% INPUTS:
%   theta: turning angle (rad)
%   M1: upstream Mach
%   gamma: ratio of specific heats
%   
% OUTPUTS:
%   M2: total downstream Mach number
%   Mn2: normal component of the Mach number after the shock
%   P0_rat: P02/P01
%   T_rat: T2/T1
%   P_rat: P2/P1
%   rho_rat: rho2/rho1
%
% Sam Jaeger
% 2/20/2024

function [M2,Mn2,P0_rat,T_rat,P_rat,rho_rat] = oblique_shock(theta,M1,gamma)
    [beta_weak, ~, ~, ~] = theta_beta_M(M1,theta,gamma);
    Mn1 = M1*sin(beta_weak);
    [Mn2,P0_rat,T_rat,P_rat,rho_rat] = normal_shock(Mn1,gamma);
    M2 = Mn2 / sin(beta_weak - theta);
end