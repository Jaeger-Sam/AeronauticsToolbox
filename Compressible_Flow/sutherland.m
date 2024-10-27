% sutherland.fcn computes the dynamic viscosity (mu) for air given a 
% temperature and using sutherlands law. If the input tmperature is above
% 1000K then will use a modified law to compute viscosity.
% https://www.cfd-online.com/Wiki/Sutherland%27s_law
% https://www.engineeringtoolbox.com/dynamic-viscosity-d_571.html
%
% INPUTS:
%   T: temperature (K or R)
%   units: 'SI' or 'US'
%
% OUTPUTS
%   mu: dynamic viscosity  (pa-s or lb/ft*s)
%
% Sam Jaeger
% 3/5/2024
%   Revised: 3/18/2024

function mu = sutherland(T,units)
    if units == 'US'
        T = T*5/9;
    end

    if T <= 0
        error('Negative temperature cannot exist!')
    elseif T<1000
        c1 = 1.458e-6;
        S = 110.4;
    
        mu = c1*((T).^(3/2))./(S + T);
    else
        T_ref = 1000; % cross over temperature
        c1 = 1.458e-6;
        S = 110.4;
        mu_ref = c1*((T_ref).^(3/2))./(S + T_ref);

        n = 0.72;
        mu = mu_ref*(T/T_ref).^n;
    end

    if units == 'US'
        mu = mu*1.49;
    end
end