% lift_coef.fcn computes the linear lift coefficient. 
%   C_L = C_L_alpha.*(alpha - alpha_L0);
%
% INPUTS:
%   alpha: angle of attack (rad)
%   C_L_alpha: lift slope (1/rad)
%   alpha_L0: (rad)
%
% OUTPUTS:
%   C_L: lift coefficient
%
% Sam Jaeger
% 1/15/2024

function C_L = lift_coef(alpha, C_L_alpha, alpha_L0)
    C_L = C_L_alpha.*(alpha - alpha_L0);
end