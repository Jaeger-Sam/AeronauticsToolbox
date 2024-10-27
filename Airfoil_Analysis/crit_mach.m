% crit_mach.fcn computes the critical Mach number given a set of
% incompressible pressure coefficients of an airfoil. User can specify
% gamma (ratio of specific heats) and the compressibility correction rule.
% See section 11.6 in Anderson Fund. of Aero. (5th Ed.). This function
% calls comp_corr.fcn as well as crit_pressure.fcn to solve the nonlinear
% system of equations.
%
% INPUTS:
%   C_p_0: vector of incompressible pressure coefficients
%   comp_corr_rule:
%       ==1 Prandtl-Gaulert
%       ==2 Karman-Tsien
%       ==3 Laitone Rule
%   gam: ratio of specific heats
%
% OUTPUTS:
%   M_cr: critical Mach number
%
% Sam Jaeger
% 1/12/2024

function M_cr = crit_mach(C_p_0,comp_corr_rule,gam)
    C_p_0_min = min(C_p_0);
    options = optimoptions('fsolve','Display','off');
    M_cr = fsolve(@(M)obj_fcn(M, C_p_0_min, comp_corr_rule, gam),0.5,options);

    function C_p_err = obj_fcn(M, C_p_0_min, comp_corr_rule, gam)
        C_p_min = comp_corr(C_p_0_min, M, comp_corr_rule,gam);
        C_p_cr = crit_pressure(M, gam);
        C_p_err = C_p_min - C_p_cr;
    end
end