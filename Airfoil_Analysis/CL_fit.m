% CL_fit.fcn computes the lift slope, angle of zero lift, maximum lift
% coefficient, and angle at maximum lift.
%   C_L = C_L_alpha.*(alpha - alpha_L0)
%
% INPUTS: 
%   alpha: vector of angles of attack (rad)
%   C_L: vector of lift coefficients
%   alpha_fit_start: angle to start of linear fit (rad)
%   alpha_fit_end: angle to end linear fit (rad)
%
% OUTPUTS:
%   C_L_alpha: lift slope (1/rad)
%   alpha_L0: angle zero lift (rad)
%   C_L_max: maximum lift coefficient
%   alpha_max: angle of attack at maximum lift coefficient
%
% Sam Jaeger
% 1/12/2024

function [C_L_alpha, alpha_L0, C_L_max, alpha_max] = CL_fit(alpha, C_L, alpha_fit_start, alpha_fit_end)
    if length(C_L) ~= length(alpha)
        error('alpha and C_L need same number of data points')
    end
    for ii=1:(length(alpha)-1)
        if alpha(ii) < alpha_fit_start && alpha(ii+1) >= alpha_fit_start
            alpha_start_index = ii;
        elseif alpha(ii) < alpha_fit_end && alpha(ii+1) >= alpha_fit_end
            alpha_end_index = ii;
        end
    end


    polycoefs = polyfit(alpha(alpha_start_index:alpha_end_index), C_L(alpha_start_index:alpha_end_index), 1);
    C_L_alpha = polycoefs(1);
    alpha_L0 = -polycoefs(2)/C_L_alpha;

    [ C_L_max, max_index] = max(C_L);
    alpha_max = alpha(max_index);
end