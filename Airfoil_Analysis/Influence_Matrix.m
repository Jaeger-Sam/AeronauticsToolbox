% P = Influence_Matrix(panel_coord_1,panel_coord_2,point)
% Influence_Matrix.fcn calculates a 2x2 influence matrix of a particular
% panel on a given point. 
% 
% The influence matrix describes the velocity induced from a single vortex 
% panel on an aribirary point in the flow and is only a function of 
% geometry. This assumes a linear distribution of vorticity between 
% panel_coord_1 and panel_coord_2. Xi and Eta are local coordinates of the 
% panel of interest.
%
% The general theory behind a vortex sheet can be found in Pg. 24 and
% 25 in Mechanics of Flight 2nd Ed. (eqns 1.5.14, 1.5.15, 1.5.16). 
% Additionally, the algorithm and notation outline for the influence matrix
% used for the vortex panel method can be found in Pg. 34 and 35 of 
% Mechanics of Flight 2nd Ed. (eqns 1.6.19 through 1.6.23).
%
% Inputs:
%   panel_coord_1: 1x2 vector of xy coordinates
%   panel_coord_2: 1x2 vector of xy coordinates
%   point: 1x2 vector of a xy point of interest
%
% Outputs:
%   P: 2x2 influence matrix 
%
% Sam Jaeger
% Written: 1/12/2022
% Revised: 12/1/2022

function P = Influence_Matrix(panel_coord_1,panel_coord_2,point)
    panel_coord_1_size = size(panel_coord_1);
    panel_coord_2_size = size(panel_coord_2);
    point_size = size(point);
    
    if panel_coord_1_size(2) ~= 2 || panel_coord_2_size(2) ~= 2 || point_size(2)  ~= 2
        error('inputs must be 1x2 column vector')
    elseif panel_coord_1_size(1) ~= 1 || panel_coord_2_size(1) ~= 1 || point_size(1)  ~= 1
        error('inputs must be 1x2 column vector')
    end
        

    dx = panel_coord_2(1) - panel_coord_1(1);
    dy = panel_coord_2(2) - panel_coord_1(2);
    l_i = sqrt( dx^2 + dy^2 );

    dx_matrix = [dx,dy;-dy,dx]; %Orientation of coordinate system
    distance_vec = point - panel_coord_1;
    distance_vec = distance_vec';
    xi_eta = (1/l_i)*dx_matrix*distance_vec;
    xi = xi_eta(1);
    eta = xi_eta(2);

    phi = atan2(eta*l_i,eta^2 + xi^2  - xi*l_i);
    psi = 0.5* log((xi^2 + eta^2) / ((xi - l_i)^2 + eta^2));

    dx_matrix_1 = [dx,-dy;dy,dx];
    xi_eta_matrix(1,1) = (l_i-xi)*phi + eta*psi;
    xi_eta_matrix(1,2) = xi*phi - eta*psi;
    xi_eta_matrix(2,1) = eta*phi - (l_i-xi)*psi - l_i;
    xi_eta_matrix(2,2) = -eta*phi - xi*psi + l_i;

    % Influence Matrix
    P = (1/(2*pi*l_i^2))*dx_matrix_1*xi_eta_matrix;
end