% cp_forces.fcn calculates the 2d (section) lift and drag coefficients
% given a set of pressure coefficient values as well as associated
% coordinates. aero_panel_geom_prep.fcn must be in the same working
% directory as cp_forces in order to work.
%
% INPUTS:
%   Cp: vector of pressure coefficients for each panel
%   coords: matrix of coordinates which describes the geometry
%   AoA: Angle of attack in deg
%
% OUTPUTS:
%   C_L: section lift coefficient
%   C_D: section drag coefficient
%   C_m_LE: section pitching moment at the leading edge
%
% Sam Jaeger
% Written: 12/17/2022
% Revised: 1/12/2024
%   Added pitching moment


function [C_L, C_D, C_m_LE] = cp_forces(Cp,coords,AoA)

    alpha = AoA*pi/180; % angle of attack in rad
    axial_airfoil_vector = [1,0]; % axial vector in airfoil coordinate 
    %                               system to compute angles

    % get panel control points, lengths, and normals
    [midpoints, dl, normal, ~]  = aero_panel_geom_prep(coords,0);
    
    % initalize data 
    force_per_panel = zeros(length(coords(:,1)),1);
    c_norm=zeros(length(coords(:,1)),1); % normal force for each panel
    c_axial=zeros(length(coords(:,1)),1); % axial force for each panel
    theta = zeros(length(coords(:,1)),1); % angle between panel normal and 
    %                                   airfoil coordinate system
    C_m_LE = 0; % initalize Pitching moment

    % loop over all points to compute forces from pressures, angles, and
    % decompose resulting forces into axial and normal components
    for ii=1:(length(midpoints(:,1)))
        % find section force per panel
        force_per_panel(ii) = Cp(ii)*dl(ii);

        % find angle between panel and airfoil coordinate system via the
        % dot product ( do not need to / by lengths because unit vectors )
        theta(ii) = acos(dot(axial_airfoil_vector,normal(ii,:)));

        % logic for angles from dot product
        if midpoints(ii,2) > 0 % top surface
            c_axial(ii) = force_per_panel(ii)*cos(theta(ii));
            c_norm(ii) = -force_per_panel(ii)*sin(theta(ii));
        elseif midpoints(ii,2) < 0 % bottom surface
            c_axial(ii) = force_per_panel(ii)*cos(-theta(ii));
            c_norm(ii) = -force_per_panel(ii)*sin(-theta(ii));
        elseif midpoints(ii,2) == 0 % if at leading or trailing edge
            c_axial(ii) = force_per_panel(ii);
            c_norm(ii) = 0;
        end

        % Pitching Moment
        C_m_LE = C_m_LE + c_norm(ii)*midpoints(ii,1) + c_axial(ii)*midpoints(ii,2);
    end

    % normal and axial coefficients
    C_N = sum(c_norm);
    C_A = sum(c_axial);
    C_m_LE = -C_m_LE; % sign convention

    % sum forces and decompose resulting axial and normal forces into
    % freestream coodinate system
    C_L = C_N*cos(alpha) - C_A*sin(alpha);
    C_D = C_N*sin(alpha) + C_A*cos(alpha);
end