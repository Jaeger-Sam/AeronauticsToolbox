% panel_alpha_sweep.fcn computes the lift, drag, pitching moment and
% pressure coefficients for a range of angle of attacks for a given airfoil
% coordinates.
%
% INPUTS:
%   alpha_deg: Mx1 vector of angle of attacks in deg
%   coords: Nx2 matrix of airfoil coordinates in Selig format
%   panel_opts: panel options data structure
%       panel_opts.V_inf: freestream velocity (m/s)
%       panel_opts.Mach: Mach number
%       panel_opts.comp_corr_rule: compressibility correction rule
%           ==1 for Prandtl-Gaulert
%           ==2 for Karman-Tsien
%           ==3 for Laitone Rule
%       panel_opts.gamma: ratio of specific heats
%       panel_opts.offset: percent offset 
%       panel_opts.kj: true/false, Uses Kutta Joukouski to compute lift and
%           pitching moment
%       panel_opts.cp: true/false, Uses Cp to compute lift and pitching 
%           moment.
%       panel_opts.critical_mach: true/false, Will compute critical mach
%           number at each angle of attack.
%       panel_opts.sweep_updates: true/false, Will display each angle of
%           attack completed.
%       panel_opts.alpha_plots: true/false
%       panel_opts.fit_data: true/false
%       panel_opts.alpha_fit_start: angle to start fitting (deg)
%       panel_opts.alpha_fit_end: angle to end fitting (deg)
%
% OUTPUTS:
%   data_out: structure with the following properties
%       C_L: Nx1 vector of section lift coefficients
%           Will retun Nx2 matrix if both Kutta-Joukouski and C_p
%           computation methods are selected. First column corresponds to
%           Kutta-Joukouski and second column corresponds to C_p.
%       C_D: Nx1 vector of section drag coefficients
%           Will retun Nx2 matrix if both Kutta-Joukouski and C_p
%           computation methods are selected. First column corresponds to
%           Kutta-Joukouski and second column corresponds to C_p.
%       C_m_LE: Nx1 vector of pitching moments
%           Will retun Nx2 matrix if both Kutta-Joukouski and C_p
%           computation methods are selected. First column corresponds to
%           Kutta-Joukouski and second column corresponds to C_p.
%       C_p: NxM matrix of pressure coefficients for each angle of attack
%       vel: NxM matrix of induced velocities around the airfoil
%       M_cr: Nx1 vector of critical Mach numbers
%       C_L_alpha: lift slope (1/rad)
%       alpha_L0: angle of zero lift (deg)
%       C_L_max: maximum lift coefficient
%       alpha_max: angle at maximum lift coefficient (deg)
%       C_D_0: drag coef. at zero lift
%       C_D_1: linear drag coef.
%       C_D_2: quadratic drag coef.
%       C_m_0: pitching moment at zero angle of attack
%       C_m_alpha: pitching moment slope (1/rad)
%
% Sam Jaeger
% 1/15/2024

function data_out = panel_alpha_sweep(alpha_deg, coords, panel_opts)
    chord_norm = coords(end,1);

    % Initalize variables
    C_p = zeros(length(coords(:,1))-1,length(alpha_deg));
    C_L_kj = zeros(length(alpha_deg),1);
    C_D_kj = zeros(length(alpha_deg),1);
    C_m_LE_kj = zeros(length(alpha_deg),1);
    C_m_c4_kj = zeros(length(alpha_deg),1);
    C_L_cp = zeros(length(alpha_deg),1);
    C_D_cp = zeros(length(alpha_deg),1);
    C_m_LE_cp = zeros(length(alpha_deg),1);
    C_m_c4_cp = zeros(length(alpha_deg),1);

    % alpha sweep
    for ii=1:length(alpha_deg)
        % vortex panel method
        % find vortex strengths
        [vortex,C_p(:,ii),vel(:,ii)] = Vortex_Panel(coords, ...
            panel_opts.Mach, ...
            panel_opts.V_inf, ...
            alpha_deg(ii)*pi/180, ...
            panel_opts.offset, ...
            panel_opts.comp_corr_rule, ...
            panel_opts.gamma);
        
        % calculate lift and pitching moment
        [C_L_kj(ii),C_m_LE_kj(ii)] = Coef_Vortex_Panel(vortex,coords,chord_norm,panel_opts.V_inf,(alpha_deg(ii)*pi/180));
        C_m_c4_kj(ii) = C_m_LE_kj(ii) + C_L_kj(ii)*chord_norm/4; % neglects drag
        
        % compare KJ calculated quantities to summed pressure
        [C_L_cp(ii), C_D_cp(ii), C_m_LE_cp(ii)] = cp_forces(C_p(:,ii),coords,alpha_deg(ii));
        C_m_c4_cp(ii) = C_m_LE_cp(ii) + C_L_cp(ii)*chord_norm/4; % neglects drag

        if panel_opts.sweep_updates == true
            disp(append('completed alpha = ',num2str(alpha_deg(ii)),' (deg)'))
        end
    end
    
    if panel_opts.kj == true && panel_opts.cp== false
        C_L = C_L_kj;
        C_m_c4 = C_m_c4_kj;
        C_D = C_D_kj;
    elseif panel_opts.cp == true && panel_opts.kj == false
        C_L = C_L_cp;
        C_m_c4 = C_m_c4_cp;
        C_D = C_D_cp;
    elseif panel_opts.kj == true && panel_opts.cp== true
        C_L(:,1) = C_L_kj;
        C_L(:,2) = C_L_cp;

        C_D(:,1) = C_D_kj;
        C_D(:,2) = C_D_cp;

        C_m_c4(:,1) = C_m_c4_kj;
        C_m_c4(:,2) = C_m_c4_cp;
    else
        error('Need to pick way to compute forces')
    end

    if panel_opts.critical_mach == true
        M_cr = zeros(length(alpha_deg),1);
        for ii=1:length(alpha_deg)
            M_cr(ii)=crit_mach(C_p(:,ii),panel_opts.comp_corr_rule,panel_opts.gamma);
        end
        data_out.M_cr = M_cr;
    end

    if panel_opts.fit_data == true
        [C_L_alpha, alpha_L0, C_L_max, alpha_max] = CL_fit(alpha_deg*pi/180, C_L, panel_opts.alpha_fit_start*pi/180, panel_opts.alpha_fit_end*pi/180);
        CL_fit_start = lift_coef(panel_opts.alpha_fit_start*pi/180, C_L_alpha, alpha_L0);
        CL_fit_end = lift_coef(panel_opts.alpha_fit_end*pi/180, C_L_alpha, alpha_L0);

        [C_D_0, C_D_1, C_D_2] = CD_fit(C_L, C_D_cp, CL_fit_start, CL_fit_end);

        [C_m_0, C_m_alpha] = Cm_fit(alpha_deg*pi/180,C_m_c4, panel_opts.alpha_fit_start*pi/180, panel_opts.alpha_fit_end*pi/180);
        
        data_out.C_L_alpha = C_L_alpha;
        data_out.alpha_L0 = alpha_L0*180/pi;
        data_out.C_L_max = C_L_max;
        data_out.alpha_max = alpha_max*180/pi;

        data_out.C_D_0 = C_D_0;
        data_out.C_D_1 = C_D_1;
        data_out.C_D_2 = C_D_2;

        data_out.C_m_0 = C_m_0;
        data_out.C_m_alpha = C_m_alpha;

        alpha_fit = linspace(alpha_deg(1),alpha_deg(end));
        C_L_fit = lift_coef(alpha_fit*pi/180,C_L_alpha,alpha_L0);
        C_D_fit = drag_coef(C_L_fit,C_D_0,C_D_1,C_D_2);
        C_m_fit = pitch_coef(alpha_fit*pi/180,C_m_0,C_m_alpha);
    end

    if panel_opts.alpha_plots == true
        if panel_opts.fit_data == false
            if size(C_L,2) <2 
                figure(521)
                subplot(3,1,1)
                plot(alpha_deg,C_L,'.')
                grid on
                xlabel('$\alpha$ (deg)','Interpreter','latex');
                ylabel('$\tilde{C_L}$','Interpreter','latex');
        
                subplot(3,1,2)
                plot(C_L,C_D,'.')
                grid on
                xlabel('$\tilde{C_L}$','Interpreter','latex');
                ylabel('$\tilde{C_D}$','Interpreter','latex');
        
                subplot(3,1,3)
                plot(alpha_deg,C_m_c4,'.')
                grid on
                ylim([min(C_m_c4)*1.1, 0])
                xlabel('$\alpha$ (deg)','Interpreter','latex');
                ylabel('$\tilde{C_m}_{c/4}$','Interpreter','latex');
            else
                figure(521)
                subplot(3,1,1)
                plot(alpha_deg,C_L(:,1),'.',alpha_deg,C_L(:,2),'*')
                grid on
                xlabel('$\alpha$ (deg)','Interpreter','latex');
                ylabel('$\tilde{C_L}$','Interpreter','latex');
                legend('KJ','Cp','Location','eastoutside')
        
                subplot(3,1,2)
                plot(C_L(:,1),C_D(:,1),'.',C_L(:,2),C_D(:,2),'*')
                grid on
                xlabel('$\tilde{C_L}$','Interpreter','latex');
                ylabel('$\tilde{C_D}$','Interpreter','latex');
                legend('KJ','Cp','Location','eastoutside')
    
                subplot(3,1,3)
                plot(alpha_deg,C_m_c4(:,1),'.',alpha_deg,C_m_c4(:,2),'*')
                grid on
                ylim([min(C_m_c4(:,1))*1.1, 0])
                xlabel('$\alpha$ (deg)','Interpreter','latex');
                ylabel('$\tilde{C_m}_{c/4}$','Interpreter','latex');
                legend('KJ','Cp','Location','eastoutside')
            end
        else % also plot fits
            if size(C_L,2) <2 
                figure(521)
                subplot(3,1,1)
                plot(alpha_deg,C_L,'.', alpha_fit, C_L_fit)
                grid on
                xlabel('$\alpha$ (deg)','Interpreter','latex');
                ylabel('$\tilde{C_L}$','Interpreter','latex');
                legend('Panel','Fit','Location','eastoutside')
        
                subplot(3,1,2)
                plot(C_L,C_D,'.', C_L_fit,C_D_fit)
                grid on
                xlabel('$\tilde{C_L}$','Interpreter','latex');
                ylabel('$\tilde{C_D}$','Interpreter','latex');
                legend('Panel','Fit','Location','eastoutside')

                subplot(3,1,3)
                plot(alpha_deg,C_m_c4,'.', alpha_fit,C_m_fit)
                grid on
                ylim([min(C_m_c4)*1.1, 0])
                xlabel('$\alpha$ (deg)','Interpreter','latex');
                ylabel('$\tilde{C_m}_{c/4}$','Interpreter','latex');
                legend('Panel','Fit','Location','eastoutside')
            else
                figure(521)
                subplot(3,1,1)
                plot(alpha_deg,C_L(:,1),'.',alpha_deg,C_L(:,2),'*', alpha_fit, C_L_fit)
                grid on
                xlabel('$\alpha$ (deg)','Interpreter','latex');
                ylabel('$\tilde{C_L}$','Interpreter','latex');
                legend('K-J','C_p','Fit','Location','eastoutside')
        
                subplot(3,1,2)
                plot(C_L(:,1),C_D(:,1),'.',C_L(:,2),C_D(:,2),'*', C_L_fit,C_D_fit)
                grid on
                xlabel('$\tilde{C_L}$','Interpreter','latex');
                ylabel('$\tilde{C_D}$','Interpreter','latex');
                legend('K-J','C_p','Fit','Location','eastoutside')
    
                subplot(3,1,3)
                plot(alpha_deg,C_m_c4(:,1),'.',alpha_deg,C_m_c4(:,2),'*', alpha_fit,C_m_fit)
                grid on
                ylim([min(C_m_c4(:,1))*1.1, 0])
                xlabel('$\alpha$ (deg)','Interpreter','latex');
                ylabel('$\tilde{C_m}_{c/4}$','Interpreter','latex');
                legend('K-J','Cp','Fit','Location','eastoutside')
            end
        end

        if panel_opts.critical_mach == true
            figure(111),plot(alpha_deg,M_cr,'o');
            xlabel('\alpha (deg)'); 
            ylabel('M_{cr}'); 
            grid on
        end
    end

    if panel_opts.cp_plots == true
        for ii=1:length(C_p(1,:))
            figure(222)
            plot(coords(1:50,1),-C_p(1:50,ii), coords(51:99,1),-C_p(51:end,ii))
            grid on
            xlabel('chord position')
            ylabel('C_p')
            leg{ii} = append('\alpha = ',num2str(alpha_deg(ii)));
            hold on
        end
        legend(leg)
        hold off
    end
    data_out.C_L = C_L;
    data_out.C_D = C_D;
    data_out.C_m_c4 = C_m_c4;
    data_out.C_p = C_p;
    data_out.vel = vel;
end