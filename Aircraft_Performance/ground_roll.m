% ground_roll.fcn computes the takeoff properties for a given aircraft and
% airport properties. Formulation from Phillips, Mechanics of Flight
% Section 3.10. Assumes aircraft operates at 70% of maximum lift. Includes
% McCormick formulation for ground effect. Assumes the lift off speed 
% (V_L0) is 1.1 times the stall speed
%
% INPUTS:
%   aircraft: data structure with properties...
%       aero.C_D_0: drag at zero lift
%       aero.C_D_1: linear drag coefficient
%       aero.C_D_2: parabolic drag coefficient
%       aero.C_L_max: maximum lift coefficient
%       aero.e: oswald efficiency factor
%       geom.h_w: main wing height above ground (ft)
%       geom.S_w: main wing area (ft^2)
%       geom.b_w: main wing span (ft^2)
%       geom.R_A: aspect ratio
%       mass.W: weight (lb)
%       propulsion.T_0: static thrust coefficient
%       propulsion.T_1: linear thrust coefficient (wrt velocity)
%       propulsion.T_2: quadratic thrust coefficient (wrt velocity)
%   mu_r: coefficient of rolling friction 
%       typically 0.04 for paved, 0.10 for grass
%   V_hw: headwind velocity (ft/s)
%   t_r: assumed rotation time (s)
%       Up to 4(s) for large transports or less than a second for small
%       aerobatic planes
%   altitude: field elevation in ft
%   t_sim: time to simulate ground roll. Needs to be long enough to lift
%       off
% 
% OUTPUTS:
%   s_g: ground roll (ft)
%   s_a: accleration distance (ft)
%   s_r: rotation distance (ft)
%   V_L0: Velocity at lift off (ft/s)
%   t_L0: time to lift off (s)
%   s_sim: simulated position (ft)
%   V_sim: simulated ground speed (ft/s)
%   t_out: simulated time (t)
% 
% Sam Jaeger
% 1/16/2024

function [s_g, s_a, s_r, V_L0, s_sim, V_sim] = ground_roll(aircraft, mu_r, V_hw, t_r, altitude, t_sim)

    rho = ATMOS(altitude,'US');
    V_L0 = 1.1*sqrt(2/aircraft.aero.C_L_max)*sqrt(aircraft.mass.W/aircraft.geom.S_w/rho);

    % integrate ode
    [t_out,X_out] = ode45(@(t,x)ground_roll_ode(t,x, aircraft, mu_r, V_hw, rho),[0,t_sim],[0,0]);
    s_sim = X_out(:,1);
    V_sim = X_out(:,2);

    % check if lift off speed is reached
    if V_sim(end) < V_L0
        error('Aircraft did not reach lift off speed, increase t_sim')
    end

    % find acceleration index
    for ii=2:length(t_out)
        if V_sim(ii-1) < V_L0 && V_sim(ii) > V_L0
            accel_index = ii;
        end
    end
    s_a = s_sim(accel_index);
    s_r = (V_L0 - V_hw)*t_r; % assuming negligible speed change during rotation
    s_g = s_a + s_r; 

    function x_dot = ground_roll_ode(t,x, aircraft, mu_r, V_hw, rho)
        g = 32.2; % accel do to gravity.

        s = x(1); % distance
        V = x(2); % ground speed

        V_air = V + V_hw; % airspeed
        %rho = ATMOS(altitude,'US');

        C_L = 0.7*aircraft.aero.C_L_max; % assume 70% of C_L_max

        % drag
        C_D = aircraft.aero.C_D_0 + ...
            aircraft.aero.C_D_1*C_L + ...
            ((16*aircraft.geom.h_w/aircraft.geom.b_w)^2)/(1 + (16*aircraft.geom.h_w/aircraft.geom.b_w)^2)*(C_L^2)/pi/aircraft.aero.e/aircraft.geom.R_A;

        D = 0.5*rho*(V_air^2)*aircraft.geom.S_w*C_D;

        % thrust
        T = aircraft.propulsion.T_0 + aircraft.propulsion.T_1*V_air + aircraft.propulsion.T_2*V_air^2; % assume parabolic dependence on airspeed
        
        % rolling friction
        F_r = mu_r*(aircraft.mass.W - 0.5*rho*(V_air^2)*aircraft.geom.S_w*C_L);

        x_dot(1) = V;
        x_dot(2) = (g/aircraft.mass.W)*(T - D - F_r);
        x_dot = x_dot';
    end
end