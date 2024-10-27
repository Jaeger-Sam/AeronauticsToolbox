% Integrand for alpha_L0. This is only vaild for 0 < x/c < p/10

% Inputs:
%   x: x location on chord
%   m: maximum camber in % of chord (camber_max = c*m/100)
%   p: location of maximum camber in tenths (x/c = p/10)
%
% Outputs:
%   dyc_dx: derivative of chord line with respect to x (single point)

function dyc_dtheta = NACA_4_dig_dyc_dx_1(theta,m,p)
    C_1 = (m/100)*(20/p);
    C_2 = (m/100)*(200/p^2);
    
    dyc_dtheta = (C_1 - C_2.*0.5.*(1 - cos(theta))).*(1 - cos(theta));
    
%    dyc_dx = (m/100)*((20/p) - (200/p^2)*(x/c));
end
