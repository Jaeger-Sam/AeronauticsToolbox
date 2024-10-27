% RK4 is a self written Runge-Kutta algorithm to integrate the equations of
% motion of FLIGHT_SIM_EOM.fcn. Described in Section 3.6 of Stengel Flight 
% Dynamics. See also Section 3.4 in Stevens and Lewis Flight Simulation and
% Control. 
%
% [sim_data_out] = RK4(aircraft,sim_options)
%
% INPUTS:
%   aircraft: data sturcture from INITALIZE_SIMULATION
%   sim_options: data structure from INITALIZE_SIMULATION
%
% OUTPUTS:
%   sim_data_out: data structure with the included variables
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
% 2/18/2024
%   Revised: 3/4/2024
%   Revised: 3/6/2024

function [sim_data_out] = RK4(aircraft,sim_options)
    ICS = sim_options.ICS;
    t = sim_options.t;
    if size(ICS) ~= 13
        error('ICS vector must be 13 x 1!')
    elseif size(t,1) > 1
        error('t must be a vector!')
    end
   
    n_step = length(t);
    %%%%%%%%%%%%%%%%%% INITALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yb = sim_options.yb.yb;
    C = aircraft.control.C;
    K = aircraft.control.K;

    x = NaN(13,n_step);
    x(:,1) = ICS;

    altitude = NaN(n_step,1); altitude(1) = - ICS(9);
    rho = NaN(n_step,1); 
    T = NaN(n_step,1); 
    p = NaN(n_step,1); 
    a = NaN(n_step,1);
    nu = NaN(n_step,1); 
    g = NaN(n_step,1);
    gb = NaN(3,n_step);
    [rho(1), T(1), p(1), a(1), nu(1), g(1), ~] = ATMOS_1976(altitude(1),'US');
    gb(:,1) = earth_to_body(x(10:13,1)',[0;0;g(1)]');

    x_dot = NaN(13,n_step); x_dot(:,1) = zeros(13,1);
    control_vec = NaN(4,n_step); control_vec(:,1) = zeros(4,1);

    phi_theta_psi = NaN(n_step,3);
    phi_theta_psi(1,:) = attitude_from_quats(ICS(10:end)')';
    
    Velocity = NaN(n_step,1); Velocity(1) = norm(ICS(1:3));
    
    F_b = NaN(n_step,3); 
    M_b  = NaN(n_step,3);
    [F_b(1,:), M_b(1,:)] = AeroForces(x(:,1), control_vec(:,1), aircraft);

    alpha  = NaN(n_step,1);
    beta_f = NaN(n_step,1);
    beta = NaN(n_step,1);
    alpha_dot = NaN(n_step,1);  alpha_dot(1) = 0;
    beta_f_dot = NaN(n_step,1); beta_f_dot(1) = 0;
    beta_dot = NaN(n_step,1); beta_dot(1) = 0;

    alpha(1) = atan2( x(3,1), x(1,1))*180/pi;
    beta_f(1) = atan2( x(2,1), x(1,1))*180/pi; % flank angle
    beta(1) = atan( cosd( alpha(1)).*tand( beta_f(1)))*180/pi; % true sideslip
    
    Mach = NaN(n_step,1); Mach(1) = Velocity(1)/a(1);
    Reynolds = NaN(n_step,1); Reynolds(1) = aircraft.geom.c_b_w*Velocity(1)/nu(1);

    %%%%%%%%%%%%%%%%%% TIME MARCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sim_options.time_integration == true
        tic
    end
    for kk=2:n_step 
        if sim_options.display_integration_time == true
            disp(append('t = ',num2str(t(kk)),' (s)'))
        end
        %%%%%%%%%%%%%%%%%%%%%%%% RK4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dt = t(kk) - t(kk-1);
        t_12 = t(kk-1) + dt/2; % midpoint of the time interval
       
        % dx1
        u = control_input(x(:,kk-1), yb, C, K, aircraft.control.max_deflect, t(kk-1));
        w = 0;
        dx1 = FLIGHT_SIM_EOM(t(:,kk-1), x(:,kk-1), u, w, aircraft, sim_options)*dt;

        % dx2
        u = control_input(x(:,kk-1) + dx1/2, yb, C, K, aircraft.control.max_deflect, t_12);
        w = 0;
        dx2 = FLIGHT_SIM_EOM(t_12, x(:,kk-1) + dx1/2, u, w, aircraft, sim_options)*dt;
        
        % dx3
        u = control_input(x(:,kk-1) + dx2/2, yb, C, K, aircraft.control.max_deflect, t_12);
        w = 0;
        dx3 = FLIGHT_SIM_EOM(t_12, x(:,kk-1) + dx2/2, u, w, aircraft, sim_options)*dt;
        
        % dx4
        u = control_input(x(:,kk-1) + dx3, yb, C, K, aircraft.control.max_deflect, t(kk));
        w = 0;
        dx4 = FLIGHT_SIM_EOM(t(kk), x(:,kk-1) + dx3, u, w, aircraft, sim_options)*dt;

        % Integrate -> Eq 3.6-16 Stengel
        x(:,kk) = x(:,kk-1) + (dx1 + 2*dx2 + 2*dx3 + dx4)/6;

        %%%%%%%%%%%%%%%%%%%% STORAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute quantities for storage
        altitude(kk) = -x(9,kk);
        [rho(kk), T(kk), p(kk), a(kk), nu(kk), g(kk), ~] = ATMOS_1976(altitude(kk),'US');
        control_vec(:,kk) = control_input(x(:,kk), yb, C, K, aircraft.control.max_deflect, t(kk));
        x_dot(:,kk) = FLIGHT_SIM_EOM(t(kk), x(:,kk), control_vec(:,kk), w, aircraft, sim_options);
        [F_b(kk,:), M_b(kk,:)] = AeroForces(x(:,kk), control_vec(:,kk), aircraft);

        Velocity(kk) = norm(x(1:3,kk));
        % aerodynamic angles (deg)
%         alpha(kk) = atan2( x(3,kk), x(1,kk))*180/pi;
%         beta_f(kk) = atan2( x(2,kk), x(1,kk))*180/pi; % flank angle
%         beta(kk) = atan( cosd( alpha(kk)).*tand( beta_f(kk)))*180/pi; % true sideslip
        [alpha(kk),beta(kk),~,beta_f(kk)] = body_to_stab(x(:,kk));

        % aero angle derivatives
        alpha_dot(kk) = (alpha(kk) - alpha(kk-1))/dt;
        beta_f_dot(kk) = (beta_f(kk) - beta_f(kk-1))/dt;
        beta_dot(kk) = (beta(kk) - beta(kk-1))/dt;

        % non dim numbers
        Mach(kk) = Velocity(kk)/a(kk);
        Reynolds(kk) = aircraft.geom.c_b_w*Velocity(kk)/nu(kk);
        
        % attitude
        phi_theta_psi(kk,:) = attitude_from_quats(x(10:13,kk)')';
        
        % body fixed accel due to gravity
        gb(:,kk) = earth_to_body(x(10:13,kk)',[0;0;g(kk)]');

        % total load factor experienced 
        n = norm(x_dot(1:3,kk))/aircraft.mass.W + 1;
        

        %%%%%%%%%%%%%%%%%%%% TERMINATE FLIGHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check state to terminate integration
        terminate = terminate_flight(t(kk),altitude(kk),alpha(kk),beta(kk),phi_theta_psi(kk,:),n,sim_options);
        if terminate == true
            break
        end

        %%%%%%%%%%%%%%%%%%%% PLOT INTEGRATED QUANTITIES %%%%%%%%%%%%%%%%%%%
        % plot relevant quantities realtime
        if sim_options.plot_real_time == true
            figure(6969)
            plot3(x(7,kk),x(8,kk),-x(9,kk),'.')
            grid on
            xlabel('$x_f$ (ft)','Interpreter','latex')
            ylabel('$y_f$ (ft)','Interpreter','latex')
            zlabel('$-z_f$ (ft)','Interpreter','latex')
            title('Integrated Positon Realtime')
            axis equal
            hold on
        end
    end
    if sim_options.plot_real_time == true
        hold off
    end
    if sim_options.time_integration == true
        toc
    end
    %%%%%%%%%%%%%%% ASSEMBLE RESULTS INTO DATA STRUCTURE %%%%%%%%%%%%%%%%%%
    % REMOVE NaN VALUES 
    sim_data_out.altitude = rmmissing(altitude)';
    sim_data_out.rho = rmmissing(rho)';
    sim_data_out.T = rmmissing(T)';
    sim_data_out.p = rmmissing(p)';
    sim_data_out.a = rmmissing(a)';
    sim_data_out.nu = rmmissing(nu)';
    sim_data_out.g = rmmissing(g)';
    sim_data_out.gb = gb(:,1:kk);
    sim_data_out.x_dot = x_dot(:,1:kk)';
    sim_data_out.control_vec = control_vec(:,1:kk)';
    sim_data_out.Velocity = rmmissing(Velocity)';
    sim_data_out.F_b = rmmissing(F_b)';
    sim_data_out.M_b  = rmmissing(M_b)';
    sim_data_out.alpha  = rmmissing(alpha)';
    sim_data_out.beta_f = rmmissing(beta_f)';
    sim_data_out.beta = rmmissing(beta)';
    sim_data_out.alpha_dot = rmmissing(alpha_dot)';
    sim_data_out.beta_f_dot = rmmissing(beta_f_dot)';
    sim_data_out.beta_dot = rmmissing(beta_dot)';
    sim_data_out.mach = rmmissing(Mach)';
    sim_data_out.Reynolds = rmmissing(Reynolds)';

    % state vector
    sim_data_out.t_out = t(1:kk);
    sim_data_out.X_out = x(:,1:kk);

    %       x(1) = u
    %       x(2) = v
    %       x(3) = w
    sim_data_out.uvw = [x(1,1:kk); x(2,1:kk); x(3,1:kk)]'/1.688; %kts
    %       x(4) = p
    %       x(5) = q
    %       x(6) = r
    sim_data_out.pqr = [x(4,1:kk); x(5,1:kk); x(6,1:kk)]'*180/pi; % deg/s
    %       x(7) = x_f
    %       x(8) = y_f
    %       x(9) = z_f
    sim_data_out.XYZ = [x(7,1:kk); x(8,1:kk); x(9,1:kk)]';
    %       x(10) = e_0
    %       x(11) = e_1
    %       x(12) = e_2
    %       x(13) = e_3
    sim_data_out.es_out = [x(10,1:kk); x(11,1:kk); x(12,1:kk); x(13,1:kk)]';

    % compute attitude
    %sim_data_out.phi_theta_psi = attitude_from_quats( x(10:13,1:kk)')*180/pi; % convert to deg
    sim_data_out.phi_theta_psi = phi_theta_psi(1:kk,:)*180/pi;
    sim_data_out.phi = sim_data_out.phi_theta_psi(1:end,1)';
    sim_data_out.theta = sim_data_out.phi_theta_psi(1:end,2)';
    sim_data_out.psi = sim_data_out.phi_theta_psi(1:end,3)';

    % flight path angle
    sim_data_out.gamma = sim_data_out.theta - sim_data_out.alpha;

    % load factors
    sim_data_out.nx_g = gb(1,1:kk)'./g(1:kk) + ( sim_data_out.x_dot(:,1))/aircraft.mass.W;
    sim_data_out.ny_g = gb(2,1:kk)'./g(1:kk) + ( sim_data_out.x_dot(:,2))/aircraft.mass.W;
    sim_data_out.nz_g = gb(3,1:kk)'./g(1:kk) + ( sim_data_out.x_dot(:,3))/aircraft.mass.W;

    % control vector
    sim_data_out.delta_T = sim_data_out.control_vec(:,1); 
    sim_data_out.delta_e = sim_data_out.control_vec(:,2);
    sim_data_out.delta_a = sim_data_out.control_vec(:,3);
    sim_data_out.delta_r = sim_data_out.control_vec(:,4);

    % earth information
    if sim_options.earth.ellipsoidal_earth == true
        sim_data_out.PHI_PSI_H = flat_to_ellipsodial(sim_data_out.t_out, sim_data_out.x_dot(7:9,:), sim_options.earth.PHI_PSI_H_0);
    elseif sim_options.earth.spherical_earth == true
        sim_data_out.PHI_PSI_H = flat_to_spherical(sim_data_out.t_out, sim_data_out.x_dot(7:9,:), sim_options.earth.PHI_PSI_H_0);
    else
        sim_data_out.PHI_PSI_H = sim_options.earth.PHI_PSI_H_0;
    end
end