% CD_fit.fcn fits parabolic coefficients given a range of lift coefficeints
% and drag coefficients to be fitted. Assumes the form of the equation...
%
%   C_D = C_D_0 + C_D_1*C_L + C_D_2*C_L^2
%
% INPUTS:
%   C_L: vector of lift coefficients 
%   C_D: vector of drag coefficients corresponding to angles of attack
%   CL_fit_start: lift coefficient to start fitting
%   CL_fit_end: lift coefficient to end fitting 
%
% OUTPUTS:
%   C_D_0: constant drag coef.
%   C_D_1: linear drag coef.
%   C_D_2: parobolic drag coef.
%
% Sam Jaeger
% 1/15/2024

function [C_D_0, C_D_1, C_D_2] = CD_fit(C_L, C_D, CL_fit_start, CL_fit_end)
    if length(C_L) ~= length(C_D)
        error('C_L and C_D need same number of data points')
    end
    for ii=1:(length(C_L)-1)
        if C_L(ii) < CL_fit_start && C_L(ii+1) >= CL_fit_start
            alpha_start_index = ii;
        elseif C_L(ii) < CL_fit_end && C_L(ii+1) >= CL_fit_end
            alpha_end_index = ii;
        end
    end

    polynomial_coefs = polyfit(C_L(alpha_start_index:alpha_end_index), C_D(alpha_start_index:alpha_end_index), 2);

    C_D_0 = polynomial_coefs(3);
    C_D_1 = polynomial_coefs(2);
    C_D_2 = polynomial_coefs(1);

end