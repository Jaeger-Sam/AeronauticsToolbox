% NewtAero.fcn returns the aerodynamic loads from Newtonian aerodynamic
% theory. Adapated from NAflightsim.m code. This code goes with the flight
% simulation toolbox. Currently only has capability for 1 control surface
% deflection.
%
% [F_b, M_b, FM] = NewtAero(x, control_vec, aircraft)
%
% INPUTS:
%   x: [13 x 1] state vector - only need u,v,w,p,q,r,z_f
%       x(1): u
%       x(2): v
%       x(3): w
%       x(4): p
%       x(5): q
%       x(6): r
%       x(7): x_f
%       x(8): y_f
%       x(9): z_f
%       x(10): e_0
%       x(11): e_x
%       x(12): e_y
%       x(13): e_z   
%   control_vec: [4x1] vector of control inputs
%       control_vec(1): delta_T
%       control_vec(2): delta_e
%       control_vec(3): delta_a
%       control_vec(4): delta_r
%   aicrcraft: data structure of aircraft with properties...
%       aero
%           modNA: logical to use modified Newtoian aerodynamic theory
%           plot_cp: logical to plot the computed cp distrbution
%       mass
%           xcm: center of mass wrt the orgin
%       geom
%           N: total number of panels
%           xc: Nx3 matrix of the center of each panel
%           n: Nx3 matrix of the normal vector of each panel
%           S: Nx3 matrix of the area of each panel
%           controls: data structure
%               n: matrix of normal vector of each control surface panel
%               xc: matrix of centers of each surface panel
%               S: matrix of areas of each panel
% 
% OUTPUTS:
%   F_b: vector of aerodynamic loads in body frame
%       F_b(1) = X
%       F_b(2) = Y
%       F_b(3) = Z
%   M_b: vector of aerodynamic moments in body frame
%       M_b(1) = MX
%       M_b(2) = MY
%       M_b(3) = MZ
%   FM: data structure of L, D, m 
%
% Sam Jaeger
% 1/29/2023
%   Revised: 2/23/2024
%   Revised: 3/4/2024
%       Added documentation
%   Revised: 3/6/2024
%       Stripped code down

function [F_b, M_b] = NewtAero(x, control_vec, aircraft)

    % rotate state vector to vehicle coordinate system
    uvw = aircraft.geom.geom_rot_mat*x(1:3); % uvw
    omega = aircraft.geom.geom_rot_mat*x(4:6); % omega

    altitude = -x(9);

    geom = aircraft.geom;
    control = aircraft.control;

    % Compute relevant quantities
    Vinf = norm(uvw); % Total Freestream Velocity
%     alpha = atan2(x(3),x(1)); % angle of attack (rad)
%     beta = atan2(x(2),Vinf); % angle of sideslip (rad)
    [rho,~,pinf,~,~,~,~] = ATMOS_1976(altitude,'SI');

    % Compute pressure coefficient and pressure distribution
    if ~aircraft.aero.modNA % Classical Newtonian Aero
        Cp0=2;
    else % modified Newtonian Aero
        gamma=1.4; % change to reacting gas eventually.
        Cp0 = 4/(gamma+1)*((gamma+1)^2/4/gamma)^(gamma/(gamma-1));
    end

    % Define relevant quantities
    qinf = 0.5*rho*Vinf^2;
    vhat = uvw/Vinf;
    N = geom.N;

    % Get direction of local velocity at each panel
    r = geom.xc - aircraft.mass.xcm';
    vloc = Vinf*vhat' - cross(r,repmat(omega',N,1),2);
    vloc = vloc./sqrt(sum(vloc.^2,2));

    % pressure coef.
    cp = sum((geom.n).*vloc,2);
    % remove contributions of non-impinging faces
    cp(cp>0)=0;
    cp = Cp0*cp.^2;
    Pdist = qinf*cp + pinf;

    % compute forces - vehicle without control
    Fdist = -Pdist.*geom.S.*geom.n;
    F=sum(Fdist,1)';
    F(abs(F)<1e-8)=0;

    % pitching moment - vehicle
    % compute moment about Center of Mass
    % positive: x- roll port down
    % y- yaw nose portward
    % z- pitch nose down
    Mdist = cross(geom.xc,Fdist,2);
    M0 = sum(Mdist,1)';
    M0(abs(M0)<1e-8)=0;

    %%%%%%%%%% CONTROLS %%%%%%%%%%%%%%%%%%%%%
    if aircraft.control.elevator == true || aircraft.control.aileron == true || aircraft.control.rudder == true
        % Define derived quantities and vectors
        delta = control_vec(2);
    
        % rotate_controls
        n_control = control.n;
        N_control = control.N;

        % rotation matrix about z axis
        R = [cos(delta) -sin(delta) 0; sin(delta) cos(delta) 0; 0 0 1];
        n_control_rot = n_control*R;
    
%         %vhat_control = [cos(alpha+delta); sin(alpha+delta); 0]; % direction of freestream
%         vhat_control_FDCS = stab_to_body(alpha+delta,beta,V_inf); % in flight dynamic coordinate system
%         vhat_control = aircraft.geom.geom_rot_mat*vhat_control_FDCS;
        vhat_control = R*vhat;

        % Get direction of local velocity at each panel
        r_control = control.xc - aircraft.mass.xcm';
        vloc = Vinf*vhat_control' - cross(r_control,repmat(omega',N_control,1),2);
        vloc = vloc./sqrt(sum(vloc.^2,2));

        cp_control = sum((n_control_rot).*vloc,2);
        % remove contributions of non-impinging faces
        cp_control(cp_control>0)=0;
        cp_control = Cp0*cp_control.^2;
        Pdist_control = qinf*cp_control + pinf;

         % compute forces - control
        Fdist_control = -Pdist_control.*control.S.*n_control_rot;
        F_control=sum(Fdist_control,1)';
        F_control(abs(F_control)<1e-8)=0;

        % pitching moment - control
        % compute moment about Center of Mass
        % positive: x- roll port down
        % y- yaw nose portward
        % z- pitch nose down
        Mdist_control = cross(control.xc,Fdist_control,2);
        M0_control = sum(Mdist_control,1)';
        M0_control(abs(M0_control)<1e-8)=0;
    else % no control
        F_control = [0;0;0];
        M0_control = [0;0;0];
    end

    %%%%%%%%%%%%%% SUM FORCES-MOMENTS %%%%%%%%%%%%%%%%
    
    % convert to standard flight dynamic CS
    F_b(1) = F(1) + F_control(1);
    F_b(2) = F(2) + F_control(2);
    F_b(3) = F(3) + F_control(3);
    F_b = F_b';

    %M_b = M0 + M0_control;
    % convert to standard flight dynamic CS
    % moment at orgin
    M_b0(1) = M0(1) + M0_control(1);
    M_b0(2) = M0(2) + M0_control(2);% - F_b(3)*aircraft.mass.xcm;
    M_b0(3) = M0(3) + M0_control(3);
    M_b0 = M_b0';

    % moment at center of mass
    %   Stengel Eq 2.3-2
    xcm = aircraft.geom.geom_rot_mat(1,:)*[aircraft.mass.xcm;  aircraft.mass.ycm;  aircraft.mass.zcm];
    ycm = aircraft.geom.geom_rot_mat(2,:)*[aircraft.mass.xcm;  aircraft.mass.ycm;  aircraft.mass.zcm];
    zcm = aircraft.geom.geom_rot_mat(3,:)*[aircraft.mass.xcm;  aircraft.mass.ycm;  aircraft.mass.zcm];
    M_b = [0, -zcm, ycm; zcm, 0, -xcm; -ycm, xcm, 0]*F_b + M_b0;

    %%%%%%%%% PLOT CP %%%%%%%%%%%%%%%%%%%%%%
    if aircraft.aero.plot_cp == true
        % actual
        alpha_a = atan2(x(3),x(1)); % angle of attack (rad)
        beta_a = atan2(x(2),norm(x(1:3))); % angle of sideslip (rad)
        if aircraft.control.elevator == true || aircraft.control.aileron == true || aircraft.control.rudder == true
            figure(6969)
            scatter3(geom.xc(:,1),geom.xc(:,2),geom.xc(:,3),40,cp,"filled")
            title(append('$C_p$ at $\alpha = $',num2str(alpha_a*180/pi),'$ \deg $',' and $\beta = $',num2str(beta_a*180/pi),'$ \deg $',' and $\delta_{bf} = $',num2str(delta*180/pi),'$ \deg $'),'Interpreter','latex','FontSize',12)
            subtitle('Displayed in geometry coordinate system','FontSize',10,'Interpreter','latex')
            xlabel('$x_{b_{cs}}$ (ft)','Interpreter','latex','FontSize',12)
            ylabel('$y_{b_{cs}}$ (ft)','Interpreter','latex','FontSize',12)
            zlabel('$z_{b_{cs}}$ (ft)','Interpreter','latex','FontSize',12)
            colorbar
            axis equal
            hold on
            scatter3(control.xc(:,1),control.xc(:,2),control.xc(:,3),40,cp_control,"filled")
            quiver3(xcm, ycm, zcm, ...
                F_b(1)/norm(F_b)*aircraft.geom.b_w*1.25, ...
                F_b(2)/norm(F_b)*aircraft.geom.b_w*1.25, ...
                F_b(3)/norm(F_b)*aircraft.geom.b_w*1.25,'LineWidth',0.8) % plot force vector
            quiver3(xcm, ycm, zcm, ...
                 M_b(1)/norm(M_b)*aircraft.geom.b_w*1.25, ...
                 M_b(2)/norm(M_b)*aircraft.geom.b_w*1.25, ...
                 M_b(3)/norm(M_b)*aircraft.geom.b_w*1.25,'LineWidth',0.8) %plot moment vector
            legend('$C_{p_{vehicle}}$','$C_{p_{control}}$','$\vec{F}_b$','$\vec{M}_b$','Interpreter','latex','Location','bestoutside')
            hold off
        else
            figure(6969)
            scatter3(geom.xc(:,1),geom.xc(:,2),geom.xc(:,3),40,cp,"filled")
            title(append('$C_p$ at $\alpha = $',num2str(alpha_a*180/pi),'$ \deg $',' and $\beta = $',num2str(beta_a*180/pi)),'Interpreter','latex','FontSize',12)
            subtitle('Displayed in geometry coordinate system','FontSize',10,'Interpreter','latex')
            xlabel('$x_{b_{cs}}$ (ft)','Interpreter','latex','FontSize',12)
            ylabel('$y_{b_{cs}}$ (ft)','Interpreter','latex','FontSize',12)
            zlabel('$z_{b_{cs}}$ (ft)','Interpreter','latex','FontSize',12)
            colorbar
            axis equal
            hold on
            quiver3(xcm, ycm, zcm, ...
                F_b(1)/norm(F_b)*aircraft.geom.b_w*1.25, ...
                F_b(2)/norm(F_b)*aircraft.geom.b_w*1.25, ...
                F_b(3)/norm(F_b)*aircraft.geom.b_w*1.25,'LineWidth',0.8) % plot force vector
            quiver3(xcm, ycm, zcm, ...
                 M_b(1)/norm(M_b)*aircraft.geom.b_w*1.25, ...
                 M_b(2)/norm(M_b)*aircraft.geom.b_w*1.25, ...
                 M_b(3)/norm(M_b)*aircraft.geom.b_w*1.25,'LineWidth',0.8) %plot moment vector
            legend('$C_{p_{vehicle}}$','$\vec{F}_b$','$\vec{M}_b$','Interpreter','latex','Location','bestoutside')
            hold off
        end
    end
    
%     % rotate back to flight dynamic coordinate system
%     F_b = aircraft.geom.geom_rot_mat*F_b;
%     M_b = aircraft.geom.geom_rot_mat*M_b;
end