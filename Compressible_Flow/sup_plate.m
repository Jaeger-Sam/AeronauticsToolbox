% [C_L_tilde, C_D_tilde] = sup_plate(alpha,M,gamma)
%
% sup_plate.m calculates the 2d lift and drag coefficients for a supersonic
% plate. See Example 9.13 in Anderson Fund. of Aero.
%
% INPUTS:
%   alpha: angle of attack (rad)
%   M: Freestream Mach number
%   gamma: ratio of specific heats
%
% OUTPUTS:
%   C_L_tilde: 2D lift coefficient
%   C_D_tilde: 2D drag coefficient
%
% Sam Jaeger
% 2/20/2024

function [C_L_tilde, C_D_tilde] = sup_plate(alpha,M,gamma)
    [~,~,~,~,P3_P1,~] = oblique_shock(alpha,M,gamma);
    [~, ~, P2_P1, ~] = expansion_wave(alpha,M,gamma);

    C_L_tilde = 2*cos(alpha)/(gamma*M^2)*(P3_P1 - P2_P1);
    C_D_tilde = 2*sin(alpha)/(gamma*M^2)*(P3_P1 - P2_P1);
end