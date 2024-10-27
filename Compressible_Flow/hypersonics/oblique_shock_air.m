% oblique_shock_air.fcn computes the properties across a oblique shock
% for air (double species diatomic gas) and includes vibrational temperature 
% effects as well as reaction effects. Reaction is assumed to be at
% equilibrium. This function calls mollier.m
%
% INPUT:
%   theta: deflection angle in rad
%   M1: Upstream Mach number
%   gamma: ratio of specific heats c_p/c_v
%   T1: Temperature (K)
%   rho1: density (kg/m^3)
%   theta_v1: vibrational temperature of species 1 (K)
%   theta_v2: vibrational temperature of species 2 (K)
%   MW1: molecular weight of species 1 (kg/kmol)
%   MW2: molecular weight of species 2 (kg/kmol)
%   y1: ratio of species 1 in a molar sense
%   y2: ratio of species 2 in a molar sense
%
% OUTPUT:
%   M2: Mach after shock
%   P0_rat: P02/P01
%   T_rat: T2/T1
%   P_rat: P2/P1
%   rho_rat: rho2/rho1
%
% Sam Jaeger
% 2/9/2024
%   Updated: 2/20/2024

function [beta,M2,T_rat,P_rat,rho_rat,u2] = oblique_shock_air(theta, M1, gamma, T1, rho1, theta_v1, theta_v2, MW1, MW2, y1, y2)
    R_univ = 8314; % J/kmol*K
    R1 = R_univ/MW1;
    R2 = R_univ/MW2;
    
    % mixture properties
    % specific energy
    e = y1*(5/2)*R1*T1 + y1*R1*theta_v1/(exp(theta_v1/T1)-1) + y2*(5/2)*R2*T1 + y2*R2*theta_v2/(exp(theta_v2/T1)-1);
    R = R1*y1 + R2*y2; % not sure if I can do this...

    % properties before shock...
    c1 = sqrt(gamma*R*T1);% Speed of sound before shock - assume perfect gas before
    u1 = M1*c1; % speed of gas
    h1 = e + R*T1; % enthalpy
    

    eps = 1/6; % inital guess for rho_1/rho_2
    tol=1e-7;
    %rho1 = .01; % kg/m^3 guess for rho1
    P1 = rho1*R*T1; % ideal gas
    T2 = 4000; % inital guess
    for ii=1:100
        rho2 = rho1/eps;
        beta = eps*tan(theta)/(1-eps) + theta;
        P2 = P1 + rho1*((u1*sin(beta))^2)*(1 - eps);
        h2 = h1 + 0.5*((u1*sin(beta))^2)*(1 - eps^2);
        
        if P2<0 || T2<0 || h2<0
            error('Pressure, Temperature, Enthalpy are negative!')
        end
        if beta<0 || beta>pi/2
            error('beta is negative or greater than 90 deg!')
        end

        % call mollier for inital gas state
        [Z2,~,~,~,T2,rho2_p] = mollier(P2,h2/1e6,0);
        
        %eps=rho1/rho2_p;
        eps_p = (P1/P2)*(T2/T1)*Z2;
        if abs(eps - eps_p) < tol
            break
        end
        eps = eps_p;
        if ii == 100
            error('Did not converge!')
        end
    
    
    end
    
    beta = eps*tan(theta)/(1-eps) + theta;
    
    % assume perfect gas after the shock
    c2 = sqrt(gamma*R*T2);
    u2_n = rho1/rho2_p *u1*sin(beta);
    %M2_n = u2_n/c2;
    u2_t = cos(beta)*u1;

    u2 = sqrt(u2_n^2 + u2_t^2);

    M2 = u2/c2;

    %M2_t = cos(beta)*M1
    %M2 = sqrt( M2_n^2 )
    
    T_rat = T2/T1;
    P_rat = P2/P1;
    rho_rat = rho2_p/rho1;
end