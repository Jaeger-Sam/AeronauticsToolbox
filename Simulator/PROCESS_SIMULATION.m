% PROCESS_SIMULATION.fcn formats the integrated simulation state vector
% from SIMULATE_FLIGHT.fcn into a data structure easy for plotting.
%
% INPUTS:
%   t_out: vector of time from integration
%   X_out: matrix of integrated state vectors
%   aircraft: data strcture of aircraft data
%   sim_options: data structure of simulation_options
%
% OUTPUTS:
%   sim_data_out: data structure with time history from flight simulation
%       including... 
%           t_out: vector of time (s)
%           X_out: matrix of integrated state variables
%           {u,v,w}: matrix of integrated air velocities (ft/s)
%           {p,q,r}: matrix of integrated angular velocities (rad/s)
%           {X,Y,Z}: matrix of integrated positions on earth (ft)
%           es_out: {e0,ex,ey,ez} quaternions
%           alpha: angle of attack (deg)
%           beta_f: flank angle (deg)
%           beta: side slip angle (deg)
%           {phi,theta,psi}, matrix of recovered attitude (deg)
%           phi: bank angle (deg)
%           theta: pitch attitude (deg)
%           psi: heading (deg)
%           gamma: flight path angle (deg)
%           altitude: height CG is above ground (ft)
%           Velocity: true airspeed (ft/s)
%           a: speed of sound (ft/s)
%           Mach: Mach number
%           Reynolds: Reynolds number
%           controls_deflection_deg: Matrix of deflections
%           x_dot: Matrix of state derivatives
%           nu: Kinematic viscosity.
%           delta_T: Vector of thrust deflections
%           delta_e: Vector of elevator deflections
%           delta_a: Vector of aileron deflections
%           delta_r: Vector of rudder deflections
%           F_b: Matrix of body fixed forces in x,y,z directions
%           M_b: Matrix of body fixed moments in x,y,z directions
%           nz_g: g load in body fixed z direction
%           PHI_PSI_H: Latitude, longitude, altitude, for flight path
%
% Sam Jaeger
% 1/15/2024

function sim_data_out = PROCESS_SIMULATION(t_out, X_out, aircraft, sim_options)
    %%%%%%% Process Data %%%%%%
    sim_data_out.t_out = t_out;
    sim_data_out.X_out = X_out;

    %       x(1) = u
    %       x(2) = v
    %       x(3) = w
    sim_data_out.uvw = [X_out(:,1), X_out(:,2), X_out(:,3)]/1.688; %kts
    %       x(4) = p
    %       x(5) = q
    %       x(6) = r
    sim_data_out.pqr = [X_out(:,4), X_out(:,5), X_out(:,6)]*180/pi; % deg/s
    %       x(7) = x_f
    %       x(8) = y_f
    %       x(9) = z_f
    sim_data_out.XYZ = [X_out(:,7), X_out(:,8), X_out(:,9)];
    %       x(10) = e_0
    %       x(11) = e_1
    %       x(12) = e_2
    %       x(13) = e_3
    sim_data_out.es_out = [X_out(:,10), X_out(:,11), X_out(:,12), X_out(:,13)];
    %       x(14) = delta_T
    %       x(15) = delta_e
    %       x(16) = delta_a
    %       x(17) = delta_r
    %sim_data_out.control_vec_out = [X_out(:,14), X_out(:,15), X_out(:,16), X_out(:,17)];

    % aerodynamic angles (deg)
    sim_data_out.alpha = atan2( sim_data_out.uvw(:,3), sim_data_out.uvw(:,1))*180/pi;
    sim_data_out.beta_f = atan2( sim_data_out.uvw(:,2), sim_data_out.uvw(:,1))*180/pi; % flank angle
    sim_data_out.beta = atan(cos( sim_data_out.alpha*pi/180).*tan( sim_data_out.beta_f*pi/180))*180/pi; % true sideslip
    
    % Attitude
    sim_data_out.phi_theta_psi = attitude_from_quats( sim_data_out.es_out)*180/pi; % convert to deg
    sim_data_out.phi = sim_data_out.phi_theta_psi(:,1);
    sim_data_out.theta = sim_data_out.phi_theta_psi(:,2);
    sim_data_out.psi = sim_data_out.phi_theta_psi(:,3);
    
    % flight path angle
    sim_data_out.gamma = sim_data_out.theta - sim_data_out.alpha;
    
    % Alitude
    sim_data_out.altitude = - sim_data_out.XYZ(:,3);
    
    % Recover Mach #, Controls, Loads
    n_steps = length(t_out); % change for variable time stepping
    sim_data_out.Velocity = zeros(n_steps,1);
    sim_data_out.a = zeros(n_steps,1);
    sim_data_out.mach = zeros(n_steps,1);
    sim_data_out.Reynolds = zeros(n_steps,1);
    sim_data_out.controls_deflection_deg = zeros(n_steps,1); % need to change this based on the model
    sim_data_out.x_dot = zeros(n_steps,13);
    sim_data_out.g = zeros(n_steps,1);
    sim_data_out.rho = zeros(n_steps,1);
    sim_data_out.T_atm = zeros(n_steps,1);
    sim_data_out.Pressure = zeros(n_steps,1);
    for ii=1:n_steps
        % Mach
        sim_data_out.Velocity(ii) = sqrt(dot( sim_data_out.uvw(ii,:), sim_data_out.uvw(ii,:)));
        [sim_data_out.rho(ii), sim_data_out.T_atm(ii), sim_data_out.Pressure(ii), sim_data_out.a(ii), sim_data_out.nu(ii), sim_data_out.g(ii),~] = ATMOS_1976(sim_data_out.altitude(ii),'US'); % speed of sound
        sim_data_out.mach(ii) = sim_data_out.Velocity(ii)/sim_data_out.a(ii);
        sim_data_out.Reynolds(ii) = aircraft.geom.c_b_w*sim_data_out.Velocity(ii)/sim_data_out.nu(ii);
    
        % controls (currently for x-15 model)
        u = control_input(X_out(ii,:)',sim_options.yb.yb,aircraft.control.C,aircraft.control.K,aircraft.control.max_deflect,t_out(ii));
        if aircraft.control.aileron == true && aircraft.control.elevator == true && aircraft.control.rudder == true && aircraft.propulsion.glider_prop_jet ~= 0
            sim_data_out.delta_T(ii,1) = u(1);
            sim_data_out.delta_e(ii,1) = u(2);
            sim_data_out.delta_a(ii,1) = u(3);
            sim_data_out.delta_r(ii,1) = u(4);
        elseif aircraft.control.aileron == false && aircraft.control.elevator == true && aircraft.control.rudder == false && aircraft.propulsion.glider_prop_jet ~= 0
            sim_data_out.delta_T(ii,1) = u(1);
            sim_data_out.delta_e(ii,1) = u(2);
            sim_data_out.delta_a(ii,1) = 0;
            sim_data_out.delta_r(ii,1) = 0;
        elseif aircraft.control.aileron == false && aircraft.control.elevator == true && aircraft.control.rudder == false && aircraft.propulsion.glider_prop_jet == 0
            sim_data_out.delta_T(ii,1) = 0;
            sim_data_out.delta_e(ii,1) = u(2);
            sim_data_out.delta_a(ii,1) = 0;
            sim_data_out.delta_r(ii,1) = 0;
        end

        % for X-15 model
        % [Right Aileron, Left Aileron, Right Stabilator, Left Stabilator, Right Flap, Left Flap, Rudder]
        %sim_data_out.controls_deflection_deg(ii,:) = [ u(2), -u(2), u(3), u(3), 0, 0, u(4)];
        
        % for MURI model
        sim_data_out.controls_deflection_deg(ii,:) = [ u(2)]'*180/pi;

        % loads - need to look at this equation
        [sim_data_out.F_b(ii,:), sim_data_out.M_b(ii,:)] = AeroForces( X_out(ii,:), u, aircraft);
        sim_data_out.nz_g(ii,1) = (-sim_data_out.F_b(ii,3)*cos(sim_data_out.alpha(ii)*pi/180) + sim_data_out.F_b(ii,1)*sin(sim_data_out.alpha(ii)*pi/180) )/aircraft.mass.W - cos(sim_data_out.phi(ii)*pi/180)*cos(sim_data_out.theta(ii)*pi/180);
    
        % Derivatives
        sim_data_out.x_dot(ii,:) = FLIGHT_SIM_EOM_constraints(t_out(ii), sim_data_out.X_out(ii,:)', aircraft, sim_options);

        % aero angle derivatives
        if ii==1
            sim_data_out.alpha_dot(ii) = 0; % for now just make it zero at inital time step
            sim_data_out.beta_dot(ii) = 0;
        else
            sim_data_out.alpha_dot(ii) = (sim_data_out.alpha(ii) - sim_data_out.alpha(ii-1))/(t_out(ii)-t_out(ii-1));
            sim_data_out.beta_dot(ii) = (sim_data_out.beta(ii) - sim_data_out.beta(ii-1))/(t_out(ii)-t_out(ii-1));
        end
    end

    % earth information
    if sim_options.earth.ellipsoidal_earth == true
        sim_data_out.PHI_PSI_H = flat_to_ellipsodial(t_out, sim_data_out.x_dot(:,7:9), sim_options.earth.PHI_PSI_H_0);
    elseif sim_options.earth.spherical_earth == true
        sim_data_out.PHI_PSI_H = flat_to_spherical(t_out, sim_data_out.x_dot(:,7:9), sim_options.earth.PHI_PSI_H_0);
    end
end