% [vortex,C_p,velocity] = Vortex_Panel(coords, Mach, V_inf, alpha, offset, comp_corr_rule, gamma)
%
% Vortex_Panel.fcn computes the strength of a series of vortex
% panels from a given set of points. The algorithm assumes a linear change
% in line vortex strength from one panel to another. 
% 
% General theory and notation on the Vortex panel method can be found in
% Mechanics of Flight 2nd Ed. Pg. 32-36.
% Algorithm outline can be found from Aeroacademy lectures on Vortex Panel 
% Method: https://aero-academy.org/courses/airfoil-analysis/
%
% Inputs:
%   coords: N x 2 matrix of coordinate points on the airfoil body. First column
%           is x coordinate, second is y coordnate. Must start on lower
%           surface and wrap to upper surface.
%   Mach: Mach number (must be < 0.8)
%   V_inf: Freestream velocity
%   alpha: angle of attack (rad)
%   offset: percent unit normal offset from panel to evaluate C_p
%   comp_corr_rule: compressibility correction rule
%       ==1 Prandtl-Gaulert
%       ==2 Karman-Tsien
%       ==3 Laitone Rule
%   gamma: ratio of specific heats
%   
% Outputs:
%   vortex: N-1 vector of vortex panel strength
%   C_p: N-1 vector of pressure coefficients for each panel
%   velocity: N-1 vector of velocity magnitudes for each vector
%
% Sam Jaeger
% Written: 1/12/2022
% Revised: 12/1/2022
%   Added C_p and velocity outputs & offset inputs
% Revised: 1/15/2024
%   Added compressibility corrections

function [vortex,C_p,velocity] = Vortex_Panel(coords, Mach, V_inf, alpha, offset, comp_corr_rule, gamma)
    %% Set up algorithm 
    N = length(coords);
    A = zeros(N,N); % Main Matrix
    B = zeros(N,1); % Aerodynamic vector
    
    x = coords(:,1);
    y = coords(:,2);
    
    dx = zeros(N,1);
    dy = zeros(N,1);
    
    [control_point, l, normal, ~] = aero_panel_geom_prep(coords,0);
    
    for i=1:(N-1)
        dx(i,1) = x(i+1) - x(i);
        dy(i,1) = y(i+1) - y(i);
    end
    
    %% Fill in [A] matrix
    
    for i=1:(N-1)
        for j=1:(N-1)

            %influence of panel j on control point i
            P = Influence_Matrix(coords(j,:),coords(j+1,:),control_point(i,:)); 
          
            A(i,j) = A(i,j) + ((x(i+1) - x(i))/l(i))*P(2,1) - ((y(i+1) - y(i))/l(i))*P(1,1);
            A(i,j+1) = A(i,j+1) + ((x(i+1) - x(i))/l(i))*P(2,2) - ((y(i+1) - y(i))/l(i))*P(1,2);   
            
        end
    end
    
    A(N,1) = 1;
    A(N,N) = 1;
    
    %% Fill in {B} vector
    
    for i=1:(N-1)
        B(i) = (V_inf/l(i))*( (y(i+1)-y(i))*cos(alpha) - (x(i+1) - x(i))*sin(alpha) );
    end
    
    %% Compute Vortex strength
    vortex = A\B;

    %% Compute Cp

    vortex_inf = zeros(N-1,2);
    velocity = zeros(N-1,1);
    C_p = zeros(N-1,1);
    for i=1:(N-1)
        for j=1:(N-1)
            P = Influence_Matrix(coords(j,:),coords(j+1,:), control_point(i,:) - (offset*normal(i,:)) );
            vortex_inf(j,:) = P*[vortex(j);vortex(j+1)];
        end
        total_vortex = sum(vortex_inf,1)'; % velocity inducted on ith panel due to all of the panels
        
        v = V_inf*[cos(alpha); sin(alpha)] + total_vortex;
        
        velocity(i) = sqrt( v(1)^2 + v(2)^2);
        C_p(i) = 1 - (velocity(i)/V_inf)^2;
    end

    %% compressibility correction
    C_p = comp_corr(C_p,Mach,comp_corr_rule,gamma);
    
 
end