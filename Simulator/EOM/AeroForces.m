% AeroForces.fcn return the aerodynamic forces and moments for a particular
% aircraft at a particular state. Gives in body fixed coordinates.
%
% INPUTS:
%   x: state vector
%       x(1) = u
%       x(2) = v
%       x(3) = w
%       x(4) = p
%       x(5) = q
%       x(6) = r
%       x(7) = x_f
%       x(8) = y_f
%       x(9) = z_f
%       x(10) = e_0
%       x(11) = e_x
%       x(12) = e_y
%       x(13) = e_z
%   control_vec: u input forces
%   aircraft: data structure of aircraft properties
% 
% Outputs: 
%   F_b: 3x1 vector of forces corresponding to body fixed coordinates x,y,z
%   M_b: 3x1 vector of moments corresponding to body fixed coords x,y,z
%   
% Sam Jaeger
% 9/26/2023
%   Revised 10/20/23
%   Revised 1/5/23
%   Revised 1/15/2024, added state vector
%   Revised: 3/6/2024, correct function call for NewtAero
%
%   TO DO:
%       Add capability for downwash lag terms
%       Add capbaility for exact alpha dot derivatives

function [F_b, M_b] = AeroForces(x, control_vec, aircraft)
    if aircraft.aero.Newtonian_Aerodynamics == true
        [F_b, M_b] = NewtAero(x,control_vec,aircraft);
    else
        u = x(1);
        v = x(2);
        w = x(3);
        p = x(4);
        q = x(5);
        r = x(6);
        altitude = -x(9);
    
        % Compute relevant quantities
        V = norm([u,v,w]); % Total Freestream Velocity
        alpha = atan2(w,u);
        beta = atan2(v,V);
        %[rho,~,~,a,nu] = ATMOS(altitude,'US');
        [rho,~,~,a,nu,~,~] = ATMOS_1976(altitude,'US');
        %[rho,~,~,a,nu,~,~] = ATMOS_1976(196850,'US');
        %rho = (0.000275462472119089)*0.0019403203319541; %kg/m^3 to slugs/ft^3
    
        M=V/a;% Mach Number
        Re=V*aircraft.geom.c_b_w/nu; % Reynolds Number
    
        % non dim rates
        p_b = p*aircraft.geom.b_w/2/V; % roll
        q_b = q*aircraft.geom.c_b_w/2/V; % pitch
        r_b = r*aircraft.geom.b_w/2/V; % yaw
    
        % !!!! APPROXIMATE alpha, beta derivatives !!!!!
        %   See Equation 1. 32 in Sinha and Ananthkrishnan AFDwEFC
        alpha_dot = q - ( -sin(beta)*(p*cos(alpha) + r*sin(alpha)) + q*cos(beta) );
        if alpha*180/pi == 0 || alpha*180/pi == 180 || alpha*180/pi == -180 || alpha*180/pi == 360 % singularity, hopefully doesn't go past 360
            beta_dot = - (r - (r*sin(alpha) + r*cos(alpha)))/cos(alpha);
        else
            beta_dot = (p - (cos(beta)*(p*cos(alpha) + r*sin(alpha)) + q*sin(beta)))/sin(alpha);
        end
    
    
        % Propulsive Forces
        [F_P, M_P] = PropForces(x, control_vec, aircraft.propulsion);
    
        % Aero Coefficients of Vehicle
        [C_L, C_D, C_S, C_l, C_m, C_n]= AeroCoefs(alpha, beta, alpha_dot, beta_dot, p_b, q_b, r_b, M, Re, control_vec, aircraft.aero);
        
        % (negative) axial force
        F_x_b = F_P(1) + 0.5*rho*(V^2)*aircraft.geom.S_w*(C_L*sin(alpha) - C_S*cos(alpha)*sin(beta) - C_D*cos(alpha)*cos(beta));
    
        % side force
        F_y_b = F_P(2) + 0.5*rho*(V^2)*aircraft.geom.S_w*(C_S*cos(beta) - C_D*sin(beta));
    
        % (negative) normal force
        F_z_b = F_P(3) + 0.5*rho*(V^2)*aircraft.geom.S_w*( -C_L*cos(alpha) - C_S*sin(alpha)*sin(beta) - C_D*sin(alpha)*cos(beta));
    
    
        % Rolling Moment
        M_x_b = M_P(1) + 0.5*rho*(V^2)*aircraft.geom.S_w*aircraft.geom.b_w*C_l;
    
        % Pitching Moment
        M_y_b = M_P(2) + 0.5*rho*(V^2)*aircraft.geom.S_w*aircraft.geom.c_b_w*C_m;
    
        % Yawing Moment
        M_z_b = M_P(3) + 0.5*rho*(V^2)*aircraft.geom.S_w*aircraft.geom.b_w*C_n;
    
    
        F_b = [F_x_b; F_y_b; F_z_b];
        M_b = [M_x_b; M_y_b; M_z_b];
    end
end