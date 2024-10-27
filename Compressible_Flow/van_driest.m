% q_dot_w = van_driest(V,alt,R_n,T_w)
%
% van_driest.fcn computes the heat transfer rate to a spherical nose in
% hypersonic flow. Uses perfect gas normal shock relations as well as air
% properties. Uses 1976 standard atmosphere.
%
% Inputs:
%   u_inf: freestream velocity (m/s)
%   alt: altitude (m)
%   R_n: nose radius (m)
%   T_w: wall temperature (K)
%
% Sam Jaeger
% 9/15/2024

function q_dot_w = van_driest(u_inf,alt,R_n,T_w)
    kappa = 0.763; % for a sphere
    gamma = 1.4; % ratio of specfic heats.
    R = 287; % J/kg*K
    cp = 7/2*R;
    Pr = 0.72; % Prandtl Number for air

    %freestream conditions
    [rho_inf, T_inf, p_inf, a_inf, ~, ~, ~] = ATMOS_1976(alt,'SI');
    M_inf = u_inf/a_inf;
    
    % normal shock
    [M2,~,T_rat,P_rat,rho_rat] = normal_shock(M_inf,gamma);
    
    % post shock, "edge conditions"
    T_e = T_rat*T_inf; % temperature (K)
    P_e = P_rat*p_inf; % pressure (Pa)
    rho_e = rho_rat*rho_inf; % density (kg/m^3)
    mu_e = sutherland(T_e,'SI'); % viscosity
    h_e = cp*T_e; % enthalpy 
    a_e = sqrt(gamma*R*T_e); % speed of sound at edge (m/s)
    u_e = M2*a_e; % velocity at edge (m/s)
    
    % wall conditions
    h_w = cp*T_w;
    
    % Adabatic wall;
    r = sqrt(Pr);
    h_aw = h_e + 0.5*r*u_e^2;
    
    % heat flux
    q_dot_w = kappa*(Pr^-0.6)*sqrt(rho_e*mu_e)/sqrt(R_n)*((2*(P_e - p_inf)/rho_e)^(1/4))*(h_aw - h_w);

end