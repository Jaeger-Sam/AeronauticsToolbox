% pitch_coef.fcn computes the linear pitching moment coefficient.
%   C_m = alpha.*C_m_alpha + C_m_0;
%
% INPUTS:
%   alpha: angle of attack
%   C_m_0: Pitching moment coefficient at zero alpha
%   C_m_alpha:  Derivative with respect to angle of attack 
%      (should be close to zero if c/4 used)
%
% OUTPUTS:
%   C_m: pitching moment coefficient
%
% Sam Jaeger
% 1/15/2024

function C_m = pitch_coef(alpha, C_m_0, C_m_alpha)
    C_m = alpha.*C_m_alpha + C_m_0;
end