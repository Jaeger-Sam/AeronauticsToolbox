% assumed_shapes_mats.fcn will generate a mass, stiffness, and coupling
% matrices for N_w bending mode and N_theta torsion modes. These modes are
% have assumed shape functions defined in Section 5.6 of Introduction to
% Structural Dynamics and Aeroelasticity. props needs to be a data
% structure with the structural properties of the wing.
%
% INPUTS:
%   N_w: number of bending modes
%   N_theta: number of torsion modes
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
%       rho_Ip: rho*I_p (density times moment of inertia xc) for torsion
%           natural frequency calculations
%
% OUTPUTS:
%   M: mass matrix (eqn 5.130)
%   K: stiffness matrix (eqn 5.130)
%   A: coupling matrix (eqn 5.125)
%   omega_w: vector of bending natural frequencies (3.257)
%   omega_theta: vector of torsion natural frequencies (3.172)
%
% Sam Jaeger
% 7/31/2024

function [M, K, A, omega_w, omega_theta] = assumed_shapes_mats(N_w, N_theta, props)
    % define properties for ease of use;
    m = props.m;
    l = props.l;
    a = props.a;
    e = props.e;
    b = props.b;
    r = props.r;
    EI = props.EI;
    GJ = props.GJ;
    %rho_Ip = props.rho_Ip;

    x_theta = e - a;

    % compute coupling matrix [A]
    A = zeros(N_theta,N_w);
    B = zeros(N_w);
    T = zeros(N_theta);

    omega_w = zeros(N_w,1); % bending natural frequencies
    omega_theta = zeros(N_theta,1); % torsion natural frequencies

    for jj = 1:N_w % loop over bending modes
        
        % table 3.1
        if jj == 1
            alpha_j = 1.87510/l;
        elseif jj == 2
            alpha_j = 4.69409/l;
        else
            alpha_j = (2*jj - 1)*pi/2 /l;
        end
        beta_j = (cosh(alpha_j*l) + cos(alpha_j*l))/(sinh(alpha_j*l) + sin(alpha_j*l));

        PHI_j = @(y) ( cosh(alpha_j.*y) - cos(alpha_j.*y) - beta_j*( sinh(alpha_j.*y) - sin(alpha_j.*y) ) ); % 3.258
        
        B(jj,jj) = (alpha_j*l)^4; %5.131

        omega_w(jj) = (alpha_j*l)^2 *sqrt(EI/m/l^4);

        for ii =1:N_theta % loop over torsion modes

            gamma_i = pi*(ii - 0.5)/l; % 5.120
            THETA_i = @(y) sqrt(2)*sin(gamma_i .*y); % 5.121

            integrand = @(y) THETA_i(y).*PHI_j(y);
            A(ii,jj) = integral( integrand, 0, l)/l; % 5.125

            T(ii,ii) = (gamma_i*l)^2; % 5.131

            %omega_theta(ii) = gamma_i*sqrt(GJ/rho_Ip);
            omega_theta(ii) = gamma_i*sqrt(GJ/m/(b^2)/(r^2));
        end
    end
    
    % mass matrix
    M = m*l*[eye(N_w), -b*x_theta*A'; - b*x_theta*A, eye(N_theta)*(b*r)^2];

    % stiffness matrix
    K = [EI/(l^3)*B, zeros(N_w,N_theta); zeros(N_theta,N_w), GJ/l*T];

end