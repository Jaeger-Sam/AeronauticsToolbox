% [M2, T_rat, P_rat, rho_rat] = expansion_wave(theta,M1,gamma)
% 
% expanson_wave.m calculates the expansion wave properties given a turning
% angle, upstream mach, and ratio of specific heats. Assumes a perfect gas.
% Additionally, assumes the expansion is isentropic.
% See Section 9.6 of Anderson Fund. of Aero.
%
% INPUT:
%   theta: turning angle (rad)
%   M1: upstream Mach number
%   gamma: ratio of specific heats
%
% OUTPUT:
%   M2: downstream Mach number
%   T_rat: T_2 / T_1 across expansion wave
%   P_rat: P_2 / P_1 across expansion wave
%   rho_rat: rho_2 / rho_1 across expansion wave
%
% Sam Jaeger
% 2/20/2024

function [M2, T_rat, P_rat, rho_rat] = expansion_wave(theta,M1,gamma)
    nu1 = Prandtl_Meyer(M1,gamma);
    
    nu_max = 2.27; % maximum Prandtl-Meyer Expansion
    theta_max = nu_max - nu1;
    if theta_max < theta
        error(append('Downstream Mach is infinite! THETA_max = ',num2str(theta_max*180/pi),' (DEG)'))
    end
    nu2 = theta + nu1;
    M2 = fzero(@(M) PM_error(nu2,M,gamma),[1.01*M1, 40*M1] );

    if M2 < M1
        error('You cannot defy the laws of Physics! M2 < M1! Did not converge!')
    end

    [T0_T_1,P0_P_1, rho0_rho_1] = isentropic(M1,gamma);
    [T0_T_2,P0_P_2, rho0_rho_2] = isentropic(M2,gamma);

    T_rat = T0_T_1/T0_T_2;
    P_rat = P0_P_1/P0_P_2;
    rho_rat = rho0_rho_1/rho0_rho_2;

    function PME = PM_error(nu2,M,gamma)
        PME = Prandtl_Meyer(M,gamma) - nu2;
    end
end