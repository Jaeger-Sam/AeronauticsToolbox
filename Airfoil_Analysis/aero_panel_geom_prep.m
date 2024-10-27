% [control_points, panel_length, normal, tangent]  = aero_panel_geom_prep(coords,disp_plot)
%
%  aero_panel_geom_prep.fcn calculates control points (midpoints), panel
%  lengths, unit normal and unit tangential vectors from a given set of
%  coordinate points.
%
% Inputs:
%   coords: Nx2 matrix of points defining a geometry
%   disp_plot: == true to display all of the plots showing normal and
%   tangential directions
%
% Outputs:
%   control_points: Nx2 matrix of control points of each panel
%   panel_length: Nx1 vector of lengths of each panel
%   normal: Nx2 matrix of unit normal vector of each panel
%   tangent: Nx2 matrix of unit tangential vector of each panel
%
% Sam Jaeger
% Written: 12/31/2021
% Revised: 12/1/2022
% Revised: 3/2/2023, updated plot
% Revised: 2/11/2024, added logic, made plots better

function [control_points, panel_length, normal, tangent]  = aero_panel_geom_prep(coords,disp_plot)
    N = length(coords(:,1));
    
    %% control points
    x = zeros(N-1,1);
    y = zeros(N-1,1);
    
    for i=1:(N-1)
        x(i,1) = (coords(i+1,1) - coords(i,1))/2 + coords(i,1);
        y(i,1) = (coords(i+1,2) - coords(i,2))/2 + coords(i,2);
    end
    
    control_points = [x,y];
    
    %% lengths
    panel_length = zeros(N-1,1);
    for i=1:(N-1)
        panel_length(i,1) = sqrt( (coords(i+1,1) - coords(i,1))^2 + (coords(i+1,2) - coords(i,2))^2 );
    end
    
    %% tangent
    tangent = zeros(N-1,2);
    for i=1:(N-1)
        tangent(i,1) = (coords(i+1,1) - coords(i,1));
        tangent(i,2) = (coords(i+1,2) - coords(i,2));
    end
    
    tangent(:,1)  = tangent(:,1) ./ panel_length;
    tangent(:,2)  = tangent(:,2) ./ panel_length;
    
    
    
    %% normal
    normal =  -[-tangent(:,2) tangent(:,1)];
    
    %% Plot

    if disp_plot == true
          figure; hold on
          quiver(control_points(:,1),control_points(:,2),normal(:,1),normal(:,2))
          quiver(control_points(:,1),control_points(:,2),tangent(:,1),tangent(:,2))
          patch(coords(:,1),coords(:,2),'black')
          axis equal
          grid on
          title('Aero Panel Geometry','Interpreter','latex')
          subtitle(append('$N_{panel} = ',num2str(N-1),'$'),'Interpreter','latex')
          xlabel('$x/c$','Interpreter','latex')
          ylabel('$y/c$','Interpreter','latex')
          hold off
    end

    %% logic to check geometry
    for ii=1:(N-1)
        if length(ii) < 1e-8
            warning('Panel length is zero! Do not have repeated coordinates! Normals and tangents will be fucked!')
        end
    end
    
    
end