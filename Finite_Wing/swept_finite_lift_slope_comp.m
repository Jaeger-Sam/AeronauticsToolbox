% swept_finite_lift_slope_comp.fcn computes the finite wing lift slope with
% compressible corrections per Equation 11.68 in Anderson Fund. of
% Aerodynamics (5th Ed.). Uses the Prandtl-Gaulert rule for compressibility
% corrections to Hembold Equation. See also Equation 5.83. 
% All inputs and outputs are scalars.
%
% Assumptions in this equation: 
%   - Subsonic compressible (not transonic)
%   - Low aspect ratio (<4)
%   - Same cross section throughout wing
%
% INPUTS:
%   C_tilde_L_alpha: 2D lift slope (1/rad)
%   M_inf: Freestream Mach number
%   e_s: Span efficiency factor
%   R_A: Aspect ratio
%   LAMBDA: sweep at half chord (deg)
%   
% OUTPUTS:
%   CL_alpha: 3D lift slope (1/rad)
%
% Sam Jaeger
% 1/12/2024

function CL_alpha = swept_finite_lift_slope_comp(C_tilde_L_alpha, M_inf, e_s, R_A, LAMBDA)
    if M_inf > 0.8 || M_inf < 0
        error(' Mach must be 0 =< M =< 0.8 ')
    end
    if R_A > 4 ||  R_A < 0.1
        error(' Aspect Ratio must be < 4 ')
    end
    CL_alpha = C_tilde_L_alpha*cosd(LAMBDA)./(sqrt(1 - (M_inf*cosd(LAMBDA)).^2 + (C_tilde_L_alpha*cosd(LAMBDA)./(pi.*R_A)).^2) + C_tilde_L_alpha*cosd(LAMBDA)./(pi*e_s*R_A));
end