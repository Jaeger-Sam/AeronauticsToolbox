% Integrand for alpha_L0. This is only vaild for p/10 < x/c < 1

% Inputs:
%   x: x location on chord
%   m: maximum camber in % of chord (camber_max = c*m/100)
%   p: location of maximum camber in tenths (x/c = p/10)
%
% Outputs:
%   dyc_dx: derivative of chord line with respect to x (single point)

function dyc_dtheta = NACA_4_dig_dyc_dx_4(theta,m,p)
    C_3 = m*2*p/10/(100 - 20*p + p^2);
    C_4 = m*2/(100 - 20*p + p^2);

    dyc_dtheta = (C_3 - C_4.*0.5.*(1 - cos(theta))).*(cos(2.*theta) - cos(theta));

%    dyc_dx = (m/(100-20*p+p^2))*((2*p/10) - 2*(x/c));
end