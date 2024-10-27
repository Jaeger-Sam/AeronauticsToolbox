% gen_aero_force_mats.fcn compute the generalized aerodynamic force
% matrices of assumed shape motion undergoing simple harmonic motion. This
% function accompanies assumed_shapes_mats.fcn and uses Theodorsens
% function. See Equation 5.129 in Hodges, Pierce Introduction to Structural
% Dyanmics and Aeroelasticity. Freestream density and velocity are not
% included in the matrices
%
% INPUTS:
%   k: reduced frequency
%   N_w: number of bending shapes
%   N_theta: number of torsion shapes
%   A: coupling matrix from assumed_shapes_mats.fcn
%   props: data structure with the following properties
%       m: mass
%       l: length
%       a: nondimensional percent distance from the half chord locating the
%           shear center (point "P"). 
%       e: nondimensional percent distance from the half chord locating the
%           center of mass (point "C"). 
%       b: semi chord
%       r: radius of gyration
%       EI: effective bending stiffness
%       GJ: effective torsion stiffness
%
% OUTPUTS:
%   M_aero: force matrix associated with gen coords acceleration
%   C_aero: force matrix associated with gen coords velocity
%   K_aero: force matrix associated with gen coords position
% 
% Sam Jaeger
% 7/31/2024

function [M_aero, C_aero, K_aero] = gen_aero_force_mats(k,N_w,N_theta,A,props)
    % define properties for ease of use;
    l = props.l;
    a = props.a;
    b = props.b;

    % compute theodorsen's function
    C = besselh(1,2,k) / (besselh(1,2,k) + 1i* besselh(0,2,k));

    % aero mass
    M_aero = - pi*(b^2)*l*[eye(N_w), b*a*A'; b*a*A, (b^2)*(a^2 + 1/8)*eye(N_theta)];

    % aero damping
    C_aero = - pi*b*l*[ 2*C*eye(N_w), -b*(1 + 2*(0.5-a)*C)*A'; 2*b*(0.5+a)*C*A, (b^2)*(0.5-a)*(1-2*(0.5+a)*C)*eye(N_theta) ];

    % aero stiffness
    K_aero = - pi*b*l*[zeros(N_w), -2*C*A'; zeros(N_theta,N_w), -b*(1+ 2*a)*C*eye(N_theta)];

end