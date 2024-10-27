% INPUTS:
%   alpha: angle of attack in radians
%   beta: angle of sideslip in radians
%   alpha_dot: time derivative of angle of attack in rad/s
%   beta_dot: time derative of angle of sideslip in rad/s
%   p_b: non dimensional roll rate
%   q_b: non dimensional pitch rate
%   r_b: non dimensional yaw rate
%   Mach: Mach number
%   Reynolds: Reynolds number
%   control_vec: u input forces
%   aero: data structure of coefficients
%
% OUTPUTS:
%   C_L: lift coefficient
%   C_D: Drag Coefficient
%   C_S: Side Force Coefficient
%   C_l: Rolling moment Coefficient
%   C_m: Pitching moment Coefficient
%   C_n: Yawing moment Coefficient
%
% Written By:
% Sam Jaeger
% 10/20/2023
%   Revised: 1/5/2024, Added input for aircraft data structure
%   Revised: 1/6/2024, Added input for Mach and Reynolds dependence
%   Revised: 1/15/2024, Added alpha_dot and beta_dot inputs

function [C_L, C_D, C_S, C_l, C_m, C_n] = AeroCoefs(alpha, beta, alpha_dot, beta_dot, p_b, q_b, r_b, Mach, Reynolds, control_vec, aero)
    
    % right now control for 4 variables
    delta_T = control_vec(1);
    delta_e = control_vec(2);
    delta_a = control_vec(3);
    delta_r = control_vec(4);

    if aero.aero_data_mat == true
        % some conditioning so matlab doesn't break itself
        if alpha > aero.alpha_table(end)
            alpha = aero.alpha_table(end);
        elseif alpha < aero.alpha_table(1)
            alpha = aero.alpha_table(1);
        end
        if delta_e > aero.delta_table(end)
            delta_e = aero.delta_table(end);
        elseif delta_e < aero.delta_table(1)
            delta_e = aero.delta_table(1);
        end

        C_L = interp2( aero.delta_table, aero.alpha_table, aero.CL, delta_e, alpha);
        C_D = interp2( aero.delta_table, aero.alpha_table, aero.CD, delta_e, alpha);
        C_m = interp2( aero.delta_table, aero.alpha_table, aero.Cm, delta_e, alpha);

        C_S =   aero.C_S_alpha*alpha...
              + aero.C_S_beta*beta...
              + aero.C_S_delta_e*delta_e...
              + aero.C_S_delta_a*delta_a...
              + aero.C_S_delta_r*delta_r;
    
        C_l =   aero.C_l_alpha*alpha...
              + aero.C_l_beta*beta...
              + aero.C_l_delta_e*delta_e...
              + aero.C_l_delta_a*delta_a...
              + aero.C_l_delta_r*delta_r;

        C_n =   aero.C_n_alpha*alpha...
              + aero.C_n_beta*beta...
              + aero.C_n_delta_e*delta_e...
              + aero.C_n_delta_a*delta_a...
              + aero.C_n_delta_r*delta_r;
    else
    
        % linear regime for now...
        C_L =   aero.C_L_0...
              + aero.C_L_alpha*alpha...
              + aero.C_L_delta_e*delta_e...
              + aero.C_L_delta_a*delta_a...
              + aero.C_L_delta_r*delta_r;
    
        if aero.drag_polar == true
            C_D =   aero.C_D_0...
                  + aero.C_D_1*C_L...
                  + aero.C_D_2*(C_L^2)...
                  + aero.C_D_delta_e*delta_e...
                  + aero.C_D_delta_a*delta_a...
                  + aero.C_D_delta_r*delta_r;
        else
            C_D =   aero.C_D_0...
                  + aero.C_D_alpha*alpha...
                  + aero.C_D_alpha2*(alpha^2)...
                  + aero.C_D_delta_e*delta_e...
                  + aero.C_D_delta_a*delta_a...
                  + aero.C_D_delta_r*delta_r;
        end
    
        C_S =   aero.C_S_alpha*alpha...
              + aero.C_S_beta*beta...
              + aero.C_S_delta_e*delta_e...
              + aero.C_S_delta_a*delta_a...
              + aero.C_S_delta_r*delta_r;
    
        C_l =   aero.C_l_alpha*alpha...
              + aero.C_l_beta*beta...
              + aero.C_l_delta_e*delta_e...
              + aero.C_l_delta_a*delta_a...
              + aero.C_l_delta_r*delta_r;
    
        C_m = aero.C_m_0...
              + aero.C_m_alpha*alpha...
              + aero.C_m_alpha2*(alpha^2)...
              + aero.C_m_beta*beta...
              + aero.C_m_delta_e*delta_e...
              + aero.C_m_delta_a*delta_a...
              + aero.C_m_delta_r*delta_r;
    
        C_n =   aero.C_n_alpha*alpha...
              + aero.C_n_beta*beta...
              + aero.C_n_delta_e*delta_e...
              + aero.C_n_delta_a*delta_a...
              + aero.C_n_delta_r*delta_r;
    end
end