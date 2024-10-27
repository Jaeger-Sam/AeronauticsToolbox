% nu = Prandtl_Meyer(M,gamma)
%
% Prandtl_Meyer.m calculates the Prandtl-Meyer function for a given mach
% number and gamma. Assumes a perfect gas. See Equation 9.42 in Anderson
% Fund. of Aero. 
%
% INPUTS:
%   M: Mach number
%   gamma: ratio of specific heats
%
% OUTPUTS:
%   nu: Prandtl-Meyer function
%
% Sam Jaeger
% 2/20/2024

function nu = Prandtl_Meyer(M,gamma)
    if M < 1
        error('Upstream Mach number must be supersonic!')
    elseif gamma < 1
        error('c_p/c_v must be greater than 1!')
    end

    t1 = sqrt((gamma+1)/(gamma-1));
    t2 = sqrt((gamma-1)/(gamma+1)*((M^2) - 1));
    t3 = sqrt((M^2) - 1);
    nu = t1*atan(t2) - atan(t3);
end