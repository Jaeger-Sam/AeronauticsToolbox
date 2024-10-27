% [beta_weak, beta_strong] = theta_beta_M(M1,theta,gamma)
%
% theta_beta_M computes the weak and strong shock solutions from the
% theta_beta_M oblique shock relations. See Equation 9.23 in Anderson Fund.
% of Aerodynamics. Will error out if there is no solution.
%
% INPUTS: 
%   M1: upstream mach number
%   theta: deflection angle (rad)
%   gamma: ratio of specific heats
%
% OUTPUS:
%   beta_weak: weak shock deflection angle (rad)
%   beta_strong: strong shock deflection angle (rad)
%   beta_max: maximum deflection angle for the given mach and gamma (deg)
%   t_max: maximum turning angle for the given mach and gamma (deg)
%
% Sam Jaeger
% 2/20/2024

function [beta_weak, beta_strong, beta_max, t_max] = theta_beta_M(M1,theta,gamma)
    % logic for checking inputs
    if gamma < 1 
        error('gamma cannot be less than 1!')
    elseif M1 < 1
        error('Mach number must be supersonic!')
    elseif pi/2 < delta
        error('Theta must less than 90 degrees!')
    elseif theta < 0
        error('Expansion wave! delta must be greater than zero!')
    end
    
    % check theta max
    beta_max = fminbnd(@(beta) theta_max(beta,M1,gamma),1*pi/180, 89.9*pi/180)*180/pi;
    t_max = -theta_max(beta_max*pi/180,M1,gamma)*180/pi;
    if t_max < theta*180/pi
        error(append('DETACHED SHOCK! No solution exists! THETA_MAX = ',num2str(t_max),' (DEG)'))
    elseif t_max - 1e-3 < theta*180/pi
        warning(append('MAXIMUM TURNING ANGLE! THETA_MAX = ',num2str(t_max),' (DEG)'))
        beta_weak = beta_max*pi/180;
        beta_strong = beta_max*pi/180;
        return
    end

    beta_weak = fzero(@(beta) theta_beta_M_error(beta,M1,theta,gamma) ,[.1 beta_max]*pi/180);
    beta_strong  =  fzero(@(beta) theta_beta_M_error(beta,M1,theta,gamma) ,[beta_max 89.9]*pi/180);

    function es = theta_beta_M_error(beta, M1, theta, gamma)
        es = ( 2*cot(beta)*( ( ((M1*sin(beta))^2) -1 )/((M1^2)*(gamma + cos(2*beta)) + 2 )) - tan(theta));
    end
    function t = theta_max(beta,M1,gamma)
        t =  - atan( 2*cot(beta)*( ( ((M1*sin(beta))^2) -1 )/((M1^2)*(gamma + cos(2*beta)) + 2 )) );
    end
end