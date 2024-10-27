% theodorsen.fcn computes the aerodynamic influence coefficients of a 2DOF
% airfoil given a reduced frequency. The general expressions for each of
% the outputs can be found on page 226 in Problem 3 of Chapter 5 of Hodges,
% Pierce Introduction to Structural Dynamics and Aeroelasticity. See also 
% Sections 5.5.1 and 5.6 for an overview of theodorsen's function and
% unsteady aerodynamics.
%
% Inputs:
%   k: reduced frequency
%   a: location of shear center (center of rotation) in percent semi-chord
%       When  a = -1, the shear center is at the LE. 
%       When  a = +1, the shear center is at the TE.
%       When  a =  0, the shear center is at the half chord (c/2 = b)
%           -1 <= a <= 1
%
% Outputs:
%   l_h: lift due to heave
%   l_theta: lift due to pitch
%   m_h: pitching moment due to heave
%   m_theta: pitching moment due to pitch
%
% Sam Jaeger
% 7/30/2024

function [l_h, l_theta, m_h, m_theta] = theodorsen(k,a)

    % compute theodorsen's function
    C = besselh(1,2,k) / (besselh(1,2,k) + 1i* besselh(0,2,k));

    l_h = 1 - 2*1i*C/k;
    l_theta = - a - 1i/k - 2*C/(k^2) - 2*1i*(0.5 - a)*C/k;
    m_h = - a + 2*1i*(0.5 + a)*C/k;
    m_theta = (1/8) + a^2 - 1i*(0.5 - a)/k + 2*(0.5 + a)*C/(k^2) + 2*1i*(0.25 - a^2)*C/k;
end