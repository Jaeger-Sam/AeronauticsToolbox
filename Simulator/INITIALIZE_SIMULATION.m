% INITALIZE_SIMULATION.fcn returns aircraft properties data structure and
% simulation options data structure for flight simulation given the
% aircraft properties file and simuation input file.
%
% [aircraft,sim_options] = INITIALIZE_SIMULATION(Aircraft_Input_Filename, Simulation_Input_Filename)
%
% INPUTS:
%   Aircraft_Input_Filename: string variable of file location defining the
%       aircraft being simulated.
%   Simulation_Input_Filname: string variable of file location defining the
%       simuation options.
%
% OUTPUTS:
%   aircraft: data structure containing information about the aircraft
%   sim_options: data structure containing information for the simulation
%
% Sam Jaeger
% 1/5/2024
%   Revised 1/9/2024: Added capability for non flat earth.
%   Revised 1/30/2024: Added Maziar's stl file reading capability
%   Revised 2/23/2024: Fixed stl file reading

function [aircraft,sim_options] = INITIALIZE_SIMULATION(Aircraft_Input_Filename, Simulation_Input_Filename)

    run(Aircraft_Input_Filename) % Get Aircraft Properties
    %%%%%%%% MASS PROPERTIES %%%%%%%%%
    mass.W = W;
    mass.Inertia = [I_xx, -I_xy, -I_xz; -I_xy, I_yy, -I_yz; -I_xz,-I_yz, I_zz];
    mass.inv_Inertia = inv(mass.Inertia); % Invert Inertia Matrix
    mass.gyro = [0, -h_z, h_y; h_z, 0, -h_x; -h_y, h_x, 0];
    mass.xcm = xcm; % center of mass (ft) - in flight dynamic frame
    mass.ycm = ycm; % center of mass (ft) - in flight dynamic frame
    mass.zcm = zcm; % center of mass (ft)

    aircraft.mass = mass;
    %%%%%% GEOMETRY PROPERTIES %%%%%%
    geom.S_w = S_w; % ft^2 % main wing area
    geom.b_w = b_w; % ft % wing span
    geom.c_b_w = c_b_w; % ft % mean aero chord
    geom.xnp = xnp; % neutral point (ft)
    geom.geom_rot_mat = geom_rot_mat; 

    %m2ft = 3.28084;% convert meters to ft
    if geom_from_stl_file == true
        % STL GEOMETRY PROPERTIES %
        % this is from Hemati's NAflightsim code get_geom.fcn
        % Read geometry data from file, then populate vehicle.geom accordingly
        %   .mat file must be in units of feet
        if exist([stlfile '.mat']) == 2
            N=load([stlfile '.mat'],'N');
            n=load([stlfile '.mat'],'n');
            v=load([stlfile '.mat'],'v');
            N=N.N;
            n=n.n;
            v=v.v;
        elseif exist(stlfile) == 2
            [N,v,n] = get_stl_data(stlfile,'US'); % will convert to ft
        else
            error('Geometry/STL data file does not exist as specified.')
        end
        
        % Fix coordinate System
        %   Commented code tries to rotate the actual geometry
        aircraft.geom.geom_rot_mat = geom_rot_mat;

%         geom.N = N;
%         geom.n = geom_rot_mat*n'; % normals
% 
%         v(:,1) = geom_rot_mat(1,1)*v(:,1) +  geom_rot_mat(1,2)*v(:,2) + geom_rot_mat(1,3)*v(:,3);
%         v(:,2) = geom_rot_mat(2,1)*v(:,1) +  geom_rot_mat(2,2)*v(:,2) + geom_rot_mat(2,3)*v(:,3);
%         v(:,3) = geom_rot_mat(3,1)*v(:,1) +  geom_rot_mat(3,2)*v(:,2) + geom_rot_mat(3,3)*v(:,3);
% 
%         v(:,4) = geom_rot_mat(1,1)*v(:,4) +  geom_rot_mat(1,2)*v(:,5) + geom_rot_mat(1,3)*v(:,6);
%         v(:,5) = geom_rot_mat(2,1)*v(:,4) +  geom_rot_mat(2,2)*v(:,5) + geom_rot_mat(2,3)*v(:,6);
%         v(:,6) = geom_rot_mat(3,1)*v(:,4) +  geom_rot_mat(3,2)*v(:,5) + geom_rot_mat(3,3)*v(:,6);
% 
%         v(:,7) = geom_rot_mat(1,1)*v(:,7) +  geom_rot_mat(1,2)*v(:,8) + geom_rot_mat(1,3)*v(:,9);
%         v(:,8) = geom_rot_mat(2,1)*v(:,7) +  geom_rot_mat(2,2)*v(:,8) + geom_rot_mat(2,3)*v(:,9);
%         v(:,9) = geom_rot_mat(3,1)*v(:,7) +  geom_rot_mat(3,2)*v(:,8) + geom_rot_mat(3,3)*v(:,9);
%    
% 
%         figure(5249)
%         plot3(v(:,1),v(:,2),v(:,3),'.'); 
%         axis equal; 
%         xlabel('$-x_b$ (ft)','Interpreter','latex'); 
%         ylabel('$y_b$ (ft)','Interpreter','latex');
%         zlabel('$-z_b$ (ft)','Interpreter','latex'); 
%         title('MURI HGV Geometry','Interpreter','latex')
% 
%         n = geom.n';
%         geom.n = n;
        
        % Given vertices v = [v1(1:3) v2(1:3) v3(1:3)]
        % determine element areas S and centroids vc
        vc = (v(:,1:3) + v(:,4:6) + v(:,7:9))./3;
        
        N=length(vc);
        S = NaN(N,1);
        for ii = 1:N
            r1 = v(ii,1:3) - v(ii,7:9);
            r2 = v(ii,1:3) - v(ii,4:6);

            % check if normals are pointing outward...
            r = vc(ii,:) - [xcm;0;0]';
            dot_prod = dot(r,n(ii,:));

            if dot_prod<0
                geom.n(ii,:) = -n(ii,:);
            end

            S(ii) = 0.5*norm(cross(r1,r2)); % panel area
        end
        geom.N = N;
        geom.n = n;
        geom.v = v;
        geom.xc = vc; % control points
        geom.S = S;
        geom.Swet = sum(S); % wetted area
        croot = max(v(:,1))-min(v(:,1));
        b = max(v(:,3))-min(v(:,3));
        
        %vehicle.geom.Sref = 0.5*b*croot; % reference area = projected area
        geom.S_w = 0.5*b*croot; % reference area = wetted area
        geom.cref = croot; % distance from nose to tail [ft]
        geom.b_w = b; % span (tip to tip) [ft]
        geom.c_b_w = 2*croot^2/3/b; % reference chord (mean areadynamic chord for delta wing, zero taper ratio)
        geom.c_r = croot;
    end

    aircraft.geom = geom;

    %%%%%%% AERODYNAMIC COEFFICIENTS %%%%%
    aircraft.aero.aero_data_mat = aero_data_mat;
    
    if aero_data_mat == true
        aircraft.aero.aero_data_matfile = aero_data_matfile;
        load(aero_data_matfile,'CL','CD','CM','alpha','delta')

%         aircraft.aero.CL = CL(:,:,7);
%         aircraft.aero.CD = CD(:,:,7);
%         aircraft.aero.Cm = CM(:,:,7);
        aircraft.aero.CL = CL(:,:);
        aircraft.aero.CD = CD(:,:);
        aircraft.aero.Cm = CM(:,:);

        aircraft.aero.alpha_table = alpha*pi/180; % be careful about rad/deg
        aircraft.aero.delta_table = delta*pi/180; 
    end

    % Specifiy aerodynamic model
    aircraft.aero.plot_cp = plot_cp;
    aircraft.aero.Newtonian_Aerodynamics = Newtonian_Aerodynamics;
    aircraft.aero.modNA = modNA;

    % lift 
    aircraft.aero.C_L_0 = C_L_0;
    aircraft.aero.C_L_alpha = C_L_alpha;

    % drag 
    aircraft.aero.drag_polar = drag_polar;
    aircraft.aero.C_D_0 = C_D_0;
    aircraft.aero.C_D_1 = C_D_1;
    aircraft.aero.C_D_2 = C_D_2;
    aircraft.aero.C_D_alpha = C_D_alpha;
    aircraft.aero.C_D_alpha2 = C_D_alpha2;
    
    % side force
    aircraft.aero.C_S_alpha = C_S_alpha;
    aircraft.aero.C_S_beta = C_S_beta;
    
    % rolling moment 
    aircraft.aero.C_l_alpha = C_l_alpha;
    aircraft.aero.C_l_beta = C_l_beta;
    
    % pitching moment
    aircraft.aero.C_m_0 = C_m_0;
    aircraft.aero.C_m_alpha = C_m_alpha;
    aircraft.aero.C_m_alpha2 = C_m_alpha2;
    aircraft.aero.C_m_beta = C_m_beta;
    
    % yawing moment
    aircraft.aero.C_n_alpha = C_n_alpha;
    aircraft.aero.C_n_beta = C_n_beta;

    % elevator
    aircraft.aero.C_L_delta_e = C_L_delta_e;
    aircraft.aero.C_D_delta_e = C_D_delta_e;
    aircraft.aero.C_S_delta_e = C_S_delta_e;
    aircraft.aero.C_l_delta_e = C_l_delta_e;
    aircraft.aero.C_m_delta_e = C_m_delta_e;
    aircraft.aero.C_n_delta_e = C_n_delta_e;

    % aileron
    aircraft.aero.C_L_delta_a = C_L_delta_a;
    aircraft.aero.C_D_delta_a = C_D_delta_a;
    aircraft.aero.C_S_delta_a = C_S_delta_a;
    aircraft.aero.C_l_delta_a = C_l_delta_a;
    aircraft.aero.C_m_delta_a = C_m_delta_a;
    aircraft.aero.C_n_delta_a = C_n_delta_a;

    % rudder
    aircraft.aero.C_L_delta_r = C_L_delta_r;
    aircraft.aero.C_D_delta_r = C_D_delta_r;
    aircraft.aero.C_S_delta_r = C_S_delta_r;
    aircraft.aero.C_l_delta_r = C_l_delta_r;
    aircraft.aero.C_m_delta_r = C_m_delta_r;
    aircraft.aero.C_n_delta_r = C_n_delta_r;

    
    %%%%%%%  PROPULSION PROPERTIES %%%%%%%
    propulsion.max_thrust = max_thrust;
    propulsion.glider_prop_jet = glider_prop_jet; 
    
    aircraft.propulsion = propulsion;
    %%%%%%% CONTROL PROPERTIES %%%%%%%%%%%
    control.full_state_info = full_state_info;
    if full_state_info == true
        control.C = eye(13,13);
    else
        control.C = C;
    end
    control.K = K;
    control.elevator = elevator;
    control.aileron = aileron;
    control.rudder = rudder;
    control.max_deflect = max_deflect*pi/180; % convert to rad

    control.xcontrol = xcontrol;
    control.ycontrol = ycontrol;

    if geom_from_stl_file == true && (elevator == true || aileron == true || rudder == true)
        xc = geom.xc;
        % % carve out a control surface in steps
        xc_control = xc;
        n_control = n;
        S_control = S;
        % carve chordwise
        idx = xc_control(:,1)>=max(xc_control(:,1))-xcontrol;
        xc_control = xc_control(idx,:);
        n_control = n_control(idx,:);
        S_control = S_control(idx);
        % carve spanwise (specified via vertical threshold in vehicle)
        idx = abs(xc_control(:,2))<=ycontrol; %0.5*zcontrol;
        xc_control = xc_control(idx,:);
        n_control = n_control(idx,:);
        S_control = S_control(idx);
        
        control.xc = xc_control;
        control.n = n_control;
        control.S = S_control;
        
        % update vehicle geom accordingly (to remove redundancies with control surface elements)
        idx=find(~ismember(xc,xc_control,'rows'));
        xc = xc(idx,:);
        n = n(idx,:);
        S = S(idx,:);
        aircraft.geom.xc = xc;
        aircraft.geom.n = n;
        aircraft.geom.S = S;
        
        aircraft.geom.N = size(xc,1);
        control.N = size(xc_control,1);
        
        %plot3(xc_control(:,1),xc_control(:,2),xc_control(:,3),'b.'); hold on;
        %axis equal
        
        xhinge = xc_control(xc_control(:,1)==min(xc_control(:,1)),:);
        xhinge = xhinge(xhinge(:,2)==min(xhinge(:,2)),:);
        rhinge = xhinge(1,:)-xhinge(2,:);
        xhinge = mean(xhinge,1);
        %plot3(xhinge(:,1),xhinge(:,2),xhinge(:,3),'go')
        %plot3([xhinge(:,1) xhinge(:,1)+rhinge(:,1)],[xhinge(:,2) xhinge(:,2)+rhinge(:,2)],[xhinge(:,3) xhinge(:,3)+rhinge(:,3)],'r-o')
        
        control.xhinge = xhinge;
        control.rhinge = rhinge;
    end


    aircraft.control = control;
    
    %%%%%%%%%%%% SIMULATION INPUTS %%%%%%%%%%%%%%%%
    run(Simulation_Input_Filename)
    %%%%%%%%%%%%%%% TIME %%%%%%%%%%%%%%%%%%%%%%%%%%
    sim_options.t_start = t_start; % s
    sim_options.t_end = t_end; % s
    sim_options.variable_dt = variable_dt;
    if variable_dt == false
        sim_options.dt = dt; 
        sim_options.n_steps = round((t_end - t_start)/dt);
        sim_options.t = t_start:dt:t_end;
    end

    %%%%%%%%%%%% INITAL CONDITONS %%%%%%%%%%%%%%%%%
    xyz_dot_start = [u0;v0;w0]*1.688; % kts to ft/s
    xyz_start = [x0;y0;z0]; % ft
    rotation_start = [p0;q0;r0]*pi/180; % deg/s to rad/s 
    es0 = eul2quat([psi0*pi/180, theta0*pi/180, phi0*pi/180],"ZYX"); % convert attitude to quaterions

    sim_options.ICS = [xyz_dot_start; rotation_start; xyz_start; es0'];

    %%%%%%%%%%%%%%% DISTRUBANCES %%%%%%%%%%%%%%%%%%
    sim_options.disturb.wind_mag_kts = wind_mag_kts;
    sim_options.disturb.wind_heading_deg = wind_heading_deg;
    sim_options.disturb.thermal_speed_kts = thermal_speed_kts;
    
    sim_options.disturb.gusts = gusts;
    sim_options.disturb.gust_mag = gust_mag;
    sim_options.disturb.gust_dir = gust_dir; 

    kt_to_fts = 1.688; 
    if wind_heading_deg <= 180 && wind_heading_deg >= 0
        V_w = [cos((wind_heading_deg*pi/180))*wind_mag_kts*kt_to_fts;
               sin((wind_heading_deg*pi/180))*wind_mag_kts*kt_to_fts;
               thermal_speed_kts*kt_to_fts];
    elseif wind_heading_deg > 180 && wind_heading_deg <= 360
        V_w = [cos(wind_heading_deg*pi/180)*wind_mag_kts*kt_to_fts;
                sin(pi - wind_heading_deg*pi/180)*wind_mag_kts*kt_to_fts;
                thermal_speed_kts*kt_to_fts];
    else
        error('Wind heading needs to be between 0 and 360 deg!')
    end
    
    sim_options.disturb.V_w = V_w;

    %%%%%%%%%%%% GUIDANCE TRAJECTORY %%%%%%%%%%%%%%
    sim_options.guidance_mat_file = guidance_mat_file; % .mat file with variable yb
    sim_options.yb = open(guidance_mat_file);

    %%%%%%%%%%%%%%%%% CONSTRAINTS %%%%%%%%%%%%%%%%%
    sim_options.constrain_roll = constrain_roll;
    sim_options.constrain_bank = constrain_bank;
    sim_options.constrain_pitch = constrain_pitch;
    sim_options.constrain_elevation = constrain_elevation;
    sim_options.constrain_yaw = constrain_yaw;
    sim_options.constrain_azimuth = constrain_azimuth;
    
    sim_options.pure_rolling = pure_rolling;
    sim_options.pure_pitching = pure_pitching;
    sim_options.pure_yawing = pure_yawing;

    %%%%%%%%%%%%% EARTH OPTIONS %%%%%%%%%%%%%%%%%%%%%
    sim_options.earth.ellipsoidal_earth = ellipsoidal_earth;
    sim_options.earth.spherical_earth = spherical_earth;
    sim_options.earth.PHI_PSI_H_0 = [Lat_0*pi/180, Long_0*pi/180, -z0];

    %%%%%%%%%%%%% DEPARTURE %%%%%%%%%%%%%%%%%%%%%%%%
    departure_deg.alpha_depart_deg = alpha_depart_deg; % angle of attack
    departure_deg.beta_depart_deg = beta_depart_deg; % true sideslip
    departure_deg.theta_depart_deg = theta_depart_deg; % pitch attitude
    departure_deg.phi_depart_deg = phi_depart_deg; % bank angle

    sim_options.departure_deg = departure_deg;

    %%%%%%%%%%%%% INTEGRATION OPTIONS %%%%%%%%%%%%%%%
    sim_options.time_integration = time_integration;
    sim_options.int_options = int_options;
    sim_options.plot_real_time = plot_real_time;
    sim_options.display_integration_time = display_integration_time;
    
    %%%%%%%%%%%%% POST PROCESSING %%%%%%%%%%%%%%%%%%%

    % Plots
    sim_options.trajectory_3d_plot = trajectory_3d_plot;
    sim_options.state_vec_plot = state_vec_plot;
    sim_options.mach_time_plot = mach_time_plot;
    sim_options.Reynolds_time_plot = Reynolds_time_plot;
    sim_options.flight_path_angle_plot = flight_path_angle_plot;
    sim_options.aero_angle_plot = aero_angle_plot;
    sim_options.load_factor_plot = load_factor_plot;
    sim_options.controls_input_plot = controls_input_plot;
    sim_options.attitude_movie_plot = attitude_movie_plot;
    
    % Attitude Movie Options
    sim_options.model_info_file = model_info_file;
    sim_options.frame_sample_time = frame_sample_time;
    sim_options.speedx = speedx; 
    sim_options.isave_movie = isave_movie;
    sim_options.movie_file_name = movie_file_name;

end