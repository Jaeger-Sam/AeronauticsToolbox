% Cm_fit.fcn computes the pitching moment at zero alpha & pitching moment
% slope. Pitching moment slope should be close to zero if moment is taken
% about c/4.
%   C_m = C_m_alpha.*alpha + C_m_0
%
% INPUTS: 
%   alpha: vector of angles of attack (rad)
%   C_m: vector of lift coefficients
%   alpha_fit_start: angle to start of linear fit (rad)
%   alpha_fit_end: angle to end linear fit (rad)
%
% OUTPUTS:
%   C_m_0: pitching moment at alpha zero
%   C_m_alpha: lift slope (1/rad)
%
% Sam Jaeger
% 1/12/2024

function [C_m_0, C_m_alpha] = Cm_fit(alpha,Cm, alpha_fit_start, alpha_fit_end)
    if length(Cm) ~= length(alpha)
        error('alpha and C_L need same number of data points')
    end
    for ii=1:(length(alpha)-1)
        if alpha(ii) < alpha_fit_start && alpha(ii+1) >= alpha_fit_start
            alpha_start_index = ii;
        elseif alpha(ii) < alpha_fit_end && alpha(ii+1) >= alpha_fit_end
            alpha_end_index = ii;
        end
    end

    polycoefs = polyfit(alpha(alpha_start_index:alpha_end_index), Cm(alpha_start_index:alpha_end_index), 1);
    C_m_alpha = polycoefs(1);
    C_m_0 = polycoefs(2);
end