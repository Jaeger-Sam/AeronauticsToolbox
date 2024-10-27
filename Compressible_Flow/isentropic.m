% [T0_T,P0_P, rho0_rho] = isentropic(M,gamma)
% 
%   isentropic.m calculates the isentropic relations of stagnation
%   conditions to the current thermodynamic state. See Equations 8.40-8.43
%   in Anderson Fund. of Aero.
%   
%   INPUTS:
%       M: Mach number
%       gamma: ratio of specific heats
%   
%   OUTPUTS:
%       T0_T: T0/T
%       P0_P: P0/P
%       rho0_rho: rho0/rho
%
%   Sam Jaeger
%   2/20/2024

function [T0_T,P0_P, rho0_rho] = isentropic(M,gamma)
    T0_T = 1 + ((gamma-1)/2)*M^2;
    P0_P = (1 + ((gamma-1)/2)*M^2)^(gamma/(gamma-1));
    rho0_rho = (1 + ((gamma-1)/2)*M^2)^(1/(gamma-1));
end