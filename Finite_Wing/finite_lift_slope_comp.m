% finite_lift_slope_comp.fcn computes the finite wing lift slope with
% compressible corrections per Equation 11.66 in Anderson Fund. of
% Aerodynamics (5th Ed.). Uses the Prandtl-Gaulert rule for compressibility
% corrections. See also Equation 5.81. All inputs and outputs are scalars.
%
% Assumptions in this equation: 
%   - Subsonic compressible (not transonic)
%   - High aspect ratio (>6)
%   - Straight wing (no sweep)
%   - Same cross section throughout wing
%
% INPUTS:
%   C_tilde_L_alpha: 2D lift slope (1/rad)
%   M_inf: Freestream Mach number
%   e_s: Span efficiency factor
%   R_A: Aspect ratio
%   
% OUTPUTS:
%   CL_alpha: 3D lift slope (1/rad)
%
% Sam Jaeger
% 1/12/2024

function CL_alpha = finite_lift_slope_comp(C_tilde_L_alpha, M_inf, e_s, R_A)
    if M_inf > 0.8 || M_inf < 0
        error(' Mach must be 0 =< M =< 0.8 ')
    end
    if R_A < 4
        error(' Aspect Ratio must be high  ')
    end
    CL_alpha = C_tilde_L_alpha./(sqrt(1 - M_inf.^2) + C_tilde_L_alpha./(pi*e_s*R_A));
end