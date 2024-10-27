% comp_corr.fcn computes the compressibility correction given an
% incompressible pressure coefficient and the freestream Mach number. User
% can choose between the Prandtl-Gaulert, Karman-Tsien, and Laitone
% corrections. See Sections 11.4 and 11.5 of Anderson Fund. of Aero. (5th
% Ed.).
%
% INPUTS:
%   C_p_0: incompressible pressure coefficient
%   M_inf: Freestream Mach number
%   comp_corr_rule:
%       ==1 Prandtl-Gaulert
%       ==2 Karman-Tsien
%       ==3 Laitone Rule
%   gam: ratio of specific heats
%
% OUTPUTS:
%   C_p: Compressible pressure coefficient
%
% Sam Jaeger
% 1/12/2024

function C_p = comp_corr(C_p_0, M_inf, comp_corr_rule,gam)
    if M_inf > 0.8 || M_inf < 0
        error(' Must be subsonic flow ')
    end

    if comp_corr_rule == 1 % Prandtl-Gaulert
        C_p = C_p_0./sqrt(1 - (M_inf.^2));
    elseif comp_corr_rule == 2 % Karman-Tsien
        C_p = C_p_0./( sqrt(1 - M_inf.^2) + ((M_inf.^2)./(1+sqrt(1-(M_inf.^2)))).*C_p_0/2);
    elseif comp_corr_rule == 3 % Laitone's Rule
        C_p = C_p_0./( sqrt(1 - M_inf.^2) + (((M_inf.^2)./(1+(((gam)-1)/2)*(M_inf.^2))/(2*sqrt(1 - M_inf.^2))).*C_p_0));
    else
        error('comp_corr must == 1,2, or 3')
    end
end