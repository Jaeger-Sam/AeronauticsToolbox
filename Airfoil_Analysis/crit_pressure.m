% crit_pressure.fcn computes the critical pressure coefficient given the
% crtical Mach number and ratio of specific heats. See Equation 11.60 from
% Anderson Fund. of Aero. (5th Ed.). Equation is derived from isentropic 
% flow relations.
%
% INPUTS:
%   M_cr: critical mach number
%   gam: ratio of specific heats
%
% OUTPUTS:
%   C_p_cr: critical pressure coefficient
%
% Sam Jaeger
% 11/12/2024

function C_p_cr = crit_pressure(M_cr,gam)
    term1 = 2./(gam.*M_cr.^2);
    term2 = 1 + ((gam-1)/2).*M_cr.^2;
    term3 = 1 + (gam-1)/2;
    term4 = gam./(gam-1);
    
    C_p_cr = term1.*((term2./term3).^term4 - 1);
end