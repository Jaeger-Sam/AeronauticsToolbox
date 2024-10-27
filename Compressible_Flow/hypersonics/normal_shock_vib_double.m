% normal_shock_vib_single.fcn computes the properties across a normal shock
% for a double species diatomic gas and includes vibrational temperature 
% effects. Assumes non reacting gas.
%
% INPUT:
%   M1: Upstream Mach number
%   gamma: ratio of specific heats c_p/c_v
%   theta_v1: vibrational temperature of species 1 (K)
%   theta_v2: vibrational temperature of species 2 (K)
%   MW1: molecular weight of species 1 (kg/kmol)
%   MW2: molecular weight of species 2 (kg/kmol)
%   y1: ratio of species 1 in a molar sense
%   y2: ratio of species 2 in a molar sense
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

function [M2,T_rat,P_rat,rho_rat] = normal_shock_vib_double(M1, gamma, T1, theta_v1, theta_v2, MW1, MW2, y1, y2, fzero_or_newton)
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
    tol=1e-6;
    rho1 = 1; % kg/m^3 guess for rho1
    P1 = rho1*R*T1; % ideal gas
    T2 = 4000; % inital guess
    for ii=1:100
        rho2 = rho1/eps;
        P2 = P1 + rho1*(u1^2)*(1 - eps);
        h2 = h1 + 0.5*(u1^2)*(1 - eps^2);
    
        if fzero_or_newton == 1
            f =  @(T2_previous) (y1*((7/2)*R1*T2_previous + R1*theta_v1/(exp(theta_v1/T2_previous)-1)) + y2*((7/2)*R2*T2_previous + R2*theta_v2/(exp(theta_v2/T2_previous)-1)) - h2);
            T2 = fzero(f,[10,10000]);
        else
    
            for jj=1:100 % Newtons Method
                T2_previous=T2;
                
                f =  y1*((7/2)*R1*T2_previous + R1*theta_v1/(exp(theta_v1/T2_previous)-1)) + y2*((7/2)*R2*T2_previous + R2*theta_v2/(exp(theta_v2/T2_previous)-1)) - h2;
                df_dT = y1*((7/2)*R1 +  R1*((theta_v1/T2_previous)^2)*exp(theta_v1/T2_previous)/((exp(theta_v1/T2_previous)-1)^2)) + y2*((7/2)*R2 +  R2*((theta_v2/T2_previous)^2)*exp(theta_v2/T2_previous)/((exp(theta_v2/T2_previous)-1)^2));
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