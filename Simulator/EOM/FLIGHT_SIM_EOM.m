% Gives derivative of state given the current state x, control vector u and 
% disturbance vector w.
%
% Flat Eart Quaternion state space EOM. See Phillips Mechanics of Flight
% Sections 11.11, 11.8.  There are different constraint EOMs. Everything
% must be in US units (ft-lb-s).
%
% Called 
%
% INPUTS:
%   t: time
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
%   u: control_vec
%       u(1): delta_T
%       u(2): delta_e
%       u(3): delta_a
%       u(4): delta_r
%   w: disturb_vec
%   aircraft: data structure of aircraft
%   sim_options: data structure of simulation options
%   
%
% OUTPUTS:
%   x_dot: vector of derivative of state space
%   aero_angs: vector of aerodynamic angles
%       aero_angs(1): alpha (rad)
%       aero_angs(2): beta (rad)
%       aero_angs(3): beta_f (rad)
%   Loads:
%       F_b: vector of body fixed aerodynamic forces
%       M_b: vector of body fixed aerodynamic moments
%   atmos_prop:
%       T_inf: Freestream temperature
%       P_inf: Freestream pressure
%       rho_inf: Freestream density
%       nu_inf: Freestream kinematic viscosity
%       M: Freestream Mach number
%       Re: Freestream Reynolds number
%       
% Sam Jaeger
% 2/18/2024
%   Revised: 3/6/2024

function x_dot = FLIGHT_SIM_EOM(t, x, control_vec, distrub_vec, aircraft, sim_options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define variables for ease of use
    u = x(1);
    v = x(2);
    w = x(3);
    
    p = x(4);
    q = x(5);
    r = x(6);
    %omega = [p; q; r]; % rotation rate
    
    x_f = x(7);
    y_f = x(8);
    z_f = x(9);
    
    E = [x(10);x(11);x(12);x(13)]; % vector of quaterions
    E = E/norm(E); % renomalize for drift error, should have length 1.
    e_0 = E(1);
    e_x = E(2);
    e_y = E(3);
    e_z = E(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Aircraft & Enviroment Properties %%%%
    mass = aircraft.mass; % Aircraft Mass
    g = grav(-z_f); % Acceleration due to gravity
    g0 = 32.17404855643; % at sea level ft/s^2
    V_w = sim_options.disturb.V_w; % Wind in earth fixed coords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Control Input %%%%%
%    control_vec = control_input(x, sim_options.yb.yb, aircraft.control.C, aircraft.control.K,aircraft.control.max_deflect,t);
%    control_vec = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Forces from Aerodynamics & Propulsion %%%%%
    [F_b, M_b] = AeroForces(x, control_vec, aircraft);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% !Derivatives! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sim_options.constrain_roll == true
        p=0;
        omega = [p; q; r]; % rotation rate

        % Eq 11.11.1
        uvw_dot = g0/mass.W*F_b + g*[2*((e_x*e_z) - (e_y*e_0)); 2*((e_y*e_z) + (e_x*e_0)); (e_z^2 + e_0^2 - e_x^2 - e_y^2)] + [((r*v)-(q*w)); (-(r*u)); ((q*u))];
        
        % Euler's full equations {omega} x [I]*{omega}
        euler_terms(1) = 0;
        euler_terms(2) = (mass.Inertia(3,3) - mass.Inertia(1,1))*p*r + mass.Inertia(1,3)*(r^2 - p^2) + mass.Inertia(1,2)*q*r - mass.Inertia(2,3)*p*q;
        euler_terms(3) = (mass.Inertia(1,1) - mass.Inertia(2,2))*p*q + mass.Inertia(1,2)*(p^2 - q^2) + mass.Inertia(2,3)*p*r - mass.Inertia(1,3)*q*r;
        
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(1) = 0;

        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2 - e_y^2 - e_z^2), 2*(e_x*e_y - e_z*e_0),           2*(e_x*e_z + e_y*e_0);
                    2*(e_x*e_y - e_z*e_0),           (e_y^2 + e_0^2 - e_x^2 - e_z^2), 2*(e_y*e_z - e_x*e_0);
                    2*(e_x*e_z - e_y*e_0),            2*(e_y*e_z + e_x*e_0),        (e_z^2 + e_0^2 - e_x^2 - e_y^2)];
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        
        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed

        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, -e_y, -e_z; e_0, -e_z, e_y; e_z, e_0, -e_x; -e_y, e_x, e_0]*omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.constrain_bank == true
        % Don't know how to do this yet, Lagrange mutlipliers?
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.constrain_pitch == true
        q=0;
        omega = [p; q; r]; % rotation rate

        % Eq 11.11.1
        uvw_dot = g0/mass.W*F_b + g*[2*((e_x*e_z) - (e_y*e_0)); 2*((e_y*e_z) + (e_x*e_0)); (e_z^2 + e_0^2 - e_x^2 - e_y^2)] + [((r*v)); ((p*w)-(r*u)); (-(p*v))];

        % Euler's full equations {omega} x [I]*{omega}
        euler_terms(1) = mass.Inertia(2,3)*( - r^2) - mass.Inertia(1,2)*p*r;
        euler_terms(2) = 0;
        euler_terms(3) = mass.Inertia(1,2)*(p^2) + mass.Inertia(2,3)*p*r;
        
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(2) = 0;

        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2 - e_y^2 - e_z^2), 2*(e_x*e_y - e_z*e_0),           2*(e_x*e_z + e_y*e_0);
                    2*(e_x*e_y - e_z*e_0),           (e_y^2 + e_0^2 - e_x^2 - e_z^2), 2*(e_y*e_z - e_x*e_0);
                    2*(e_x*e_z - e_y*e_0),            2*(e_y*e_z + e_x*e_0),        (e_z^2 + e_0^2 - e_x^2 - e_y^2)];
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        
        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed

        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, -e_y, -e_z; e_0, -e_z, e_y; e_z, e_0, -e_x; -e_y, e_x, e_0]*omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.constrain_elevation == true
        % Don't know how to do this yet, Lagrange mutlipliers?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.constrain_yaw == true
        r=0;
        omega = [p; q; r]; % rotation rate

        % Eq 11.11.1
        uvw_dot = g0/mass.W*F_b + g*[2*((e_x*e_z) - (e_y*e_0)); 2*((e_y*e_z) + (e_x*e_0)); (e_z^2 + e_0^2 - e_x^2 - e_y^2)] + [(-(q*w)); ((p*w)); ((q*u)-(p*v))];

        % Euler's full equations {omega} x [I]*{omega}
        euler_terms(1) = mass.Inertia(2,3)*(q^2) + mass.Inertia(1,3)*p*q;
        euler_terms(2) = mass.Inertia(1,3)*(-p^2) - mass.Inertia(2,3)*p*q;
        euler_terms(3) = 0;
        
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(3) = 0;

        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2 - e_y^2 - e_z^2), 2*(e_x*e_y - e_z*e_0),           2*(e_x*e_z + e_y*e_0);
                    2*(e_x*e_y - e_z*e_0),           (e_y^2 + e_0^2 - e_x^2 - e_z^2), 2*(e_y*e_z - e_x*e_0);
                    2*(e_x*e_z - e_y*e_0),            2*(e_y*e_z + e_x*e_0),        (e_z^2 + e_0^2 - e_x^2 - e_y^2)];
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        
        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed

        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, -e_y, -e_z; e_0, -e_z, e_y; e_z, e_0, -e_x; -e_y, e_x, e_0]*omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.constrain_azimuth == true
        % Don't know how to do this yet, Lagrange mutlipliers?
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.pure_rolling == true
        q=0;
        r=0;
        omega = [p; q; r]; % rotation rate
        %e_y=0; e_z=0; E = [e_0;e_x;e_y;e_z];

        % Eq 11.11.1
        uvw_dot = g0/mass.W*F_b + g*[0; 2*( (e_x*e_0)); ( e_0^2 - e_x^2 )] + [0; ((p*w)); (-(p*v))];
        
        % Euler's equations - constraints for 2 rotations makes everything zero
        euler_terms = [0,0,0];
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(2)=0;
        omega_dot(3)=0;

        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2),              0,                0;
                    0,              (e_0^2 - e_x^2),     2*(-e_x*e_0);
                    0,                 2*( e_x*e_0), (e_0^2 - e_x^2)];
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        
        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed
        
        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, 0, 0; e_0, 0, 0; 0, e_0, -e_x; 0, e_x, e_0]*omega;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.pure_pitching == true
        p=0;
        r=0;
        omega = [p; q; r]; % rotation rate
        %e_x=0; e_z=0; E = [e_0;e_x;e_y;e_z];

        % Eq 11.11.1
        uvw_dot = g0/mass.W*F_b + g*[2*(- (e_y*e_0)); 0; ( e_0^2 - e_y^2)] + [(-(q*w)); 0; ((q*u))];
        
        % Euler's equations - constraints for 2 rotations makes everything zero
        euler_terms = [0,0,0];
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(1)=0;
        omega_dot(3)=0;

        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_0^2 - e_y^2), 0,           2*( e_y*e_0);
                    0,           (e_y^2 + e_0^2), 0;
                    2*( - e_y*e_0),            0,        ( e_0^2 - e_y^2)];
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        
        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed
        
        % Eq. 11.11.4
        E_dot = 0.5*[0, -e_y, 0; e_0, 0, e_y; 0, e_0, 0; -e_y, 0, e_0]*omega;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.pure_yawing == true
        p=0;
        q=0;
        omega = [p; q; r]; % rotation rate
        %e_x=0; e_y=0; E = [e_0;e_x;e_y;e_z];
    
        % Eq 11.11.1
        uvw_dot = g0/mass.W*F_b + g*[0; 0; (e_z^2 + e_0^2)] + [((r*v)); (-(r*u)); 0];
        
        % Euler's equations - constraints for 2 rotations makes everything zero
        euler_terms = [0,0,0];
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(1)=0;
        omega_dot(2)=0;

        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_0^2 - e_z^2), 2*( - e_z*e_0),           0;
                    2*( e_z*e_0),           ( e_0^2 - e_z^2),  0;
                    0,            0,        (e_z^2 + e_0^2)];
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        
        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed
        
        % Eq. 11.11.4
        E_dot = 0.5*[0, 0, -e_z; e_0, -e_z, 0; e_z, e_0, 0; 0, 0, e_0]*omega;
%%%%%%%%%%%%%%%%%%%%%%%%%
    else % unconstrained

        omega = [p; q; r]; % rotation rate
        % Eq 11.11.1
        uvw_dot = g0/mass.W*F_b + g*[2*((e_x*e_z) - (e_y*e_0)); 2*((e_y*e_z) + (e_x*e_0)); (e_z^2 + e_0^2 - e_x^2 - e_y^2)] + [((r*v)-(q*w)); ((p*w)-(r*u)); ((q*u)-(p*v))];
        
        % Euler's full equations {omega} x [I]*{omega}
        euler_terms(1) = (mass.Inertia(2,2) - mass.Inertia(3,3))*q*r + mass.Inertia(2,3)*(q^2 - r^2) + mass.Inertia(1,3)*p*q - mass.Inertia(1,2)*p*r;
        euler_terms(2) = (mass.Inertia(3,3) - mass.Inertia(1,1))*p*r + mass.Inertia(1,3)*(r^2 - p^2) + mass.Inertia(1,2)*q*r - mass.Inertia(2,3)*p*q;
        euler_terms(3) = (mass.Inertia(1,1) - mass.Inertia(2,2))*p*q + mass.Inertia(1,2)*(p^2 - q^2) + mass.Inertia(2,3)*p*r - mass.Inertia(1,3)*q*r;
        
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        
        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2 - e_y^2 - e_z^2), 2*(e_x*e_y - e_z*e_0),           2*(e_x*e_z + e_y*e_0);
                    2*(e_x*e_y - e_z*e_0),           (e_y^2 + e_0^2 - e_x^2 - e_z^2), 2*(e_y*e_z - e_x*e_0);
                    2*(e_x*e_z - e_y*e_0),            2*(e_y*e_z + e_x*e_0),        (e_z^2 + e_0^2 - e_x^2 - e_y^2)];
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        
        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed
        
        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, -e_y, -e_z; e_0, -e_z, e_y; e_z, e_0, -e_x; -e_y, e_x, e_0]*omega;
    end
    % State vector
    x_dot = [uvw_dot; omega_dot; xyz_dot; E_dot];
end