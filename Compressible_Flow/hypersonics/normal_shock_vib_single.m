% normal_shock_vib_single.fcn computes the properties across a normal shock
% for a single species diatomic gas and includes vibrational temperature 
% effects.
%
% INPUT:
%   M1: Upstream Mach number
%   gamma: ratio of specific heats c_p/c_v
%   theta_v: vibrational temperature (K)
%   MW: molecular weight (kg/kmol)
%   fzero_or_newton: ==1 then fzero is used to solve for T, else a newtons
%       method is used to solve for T
%
% OUTPUT:
%   M2: Mach after shock
%   P0_rat: P02/P01
%   T_rat: T2/T1
%   P_rat: P2/P1
%   rho_rat: rho2/rho1
%
% Sam Jaeger
% 2/8/2024

function [M2,T_rat,P_rat,rho_rat] = normal_shock_vib_single(M1, gamma, T1, theta_v, MW, fzero_or_newton)
    R_univ = 8314; % J/kmol*K
    R = R_univ/MW;
    
    % specific energy
    e = (5/2)*R*T1 + R*theta_v/(exp(theta_v/T1)-1);
    
    % Speed of sound before shock - assume perfect gas before
    c1 = sqrt(gamma*R*T1);
    u1 = M1*c1; % speed of gas
    h1 = e + R*T1; % enthalpy
    
    eps = 1/6; % inital guess for rho_1/rho_2
    tol=1e-6;
    rho1 = 1; % kg/m^3 guess for rho1
    P1 = rho1*R*T1; % ideal gas
    T2 = 4000; % inital guess
    for ii=1:100
        rho2 = rho1/eps;
        P2 = P1 + rho1*(u1^2)*(1 - eps);
        h2 = h1 + 0.5*(u1^2)*(1 - eps^2);
    
        if fzero_or_newton == 1
            f =  @(T2_previous) (((7/2)*R*T2_previous + R*theta_v/(exp(theta_v/T2_previous)-1)) - h2);
            T2 = fzero(f,[10,10000]);
        else
    
            for jj=1:100 % Newtons Method
                T2_previous=T2;
                
                f =  (7/2)*R*T2_previous + R*theta_v/(exp(theta_v/T2_previous)-1) - h2;
                df_dT = (7/2)*R +  R*((theta_v/T2_previous)^2)*exp(theta_v/T2_previous)/((exp(theta_v/T2_previous)-1)^2);
                T2 = T2_previous - f/df_dT;
        
                if abs(T2 - T2_previous) < tol*1e-5
                    break
                end
                if jj==100
                    error('Temperature did not converge!')
                end
            end
        end
    
        rho2_p = P2/R/T2;
        if abs(rho2 - rho2_p) < tol
            break
        end
    
        eps=rho1/rho2_p;
        if ii == 100
            error('Did not converge!')
        end
    
    
    end
    
    % assume perfect gas after the shock
    c2 = sqrt(gamma*R*T2);
    u2 = rho1/rho2_p *u1;
    M2 = u2/c2;
    
    T_rat = T2/T1;
    P_rat = P2/P1;
    rho_rat = rho2_p/rho1;
end