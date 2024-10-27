% SIMULATE_FLIGHT.fcn simulates the flight of the aircraft defined from
% aircraft and sim_options data structures. Returns all data in US units
% (ft - lb - s).
%
% INPUTS:
%   aircraft: data structure from initalize_sim.fcn
%   sim_options: data structure from initalize_sim.fcn 
%
% OUTPUTS:
%   t_out: vector of time from integration
%   X_out: matrix of integrated state vectors
%      
% Sam Jaeger
% 1/6/2024
%   Revised 1/9/2024: Added capability for non flat earth.
%   Revised 1/15/2024: Moved postprocessing to separate function

function [t_out, X_out] = SIMULATE_FLIGHT(aircraft,sim_options)

%     sim_options.int_options.Events = @(t,x)terminate_flight(t,x,sim_options.departure_deg);
    % integrate
    if sim_options.time_integration == true
        if sim_options.variable_dt == true % variable time stepping with MATLAB
            tic
            [t_out,X_out] = ode45(@(t,x)FLIGHT_SIM_EOM_constraints(t,x,aircraft,sim_options),[sim_options.t_start,sim_options.t_end],sim_options.ICS,sim_options.int_options); 
            toc
        else
            tic
            [t_out,X_out] = ode45(@(t,x)FLIGHT_SIM_EOM_constraints(t,x,aircraft,sim_options), sim_options.t, sim_options.ICS, sim_options.int_options);
            toc
        end
    else
        if sim_options.variable_dt == true % variable time stepping with MATLAB
            [t_out,X_out] = ode45(@(t,x)FLIGHT_SIM_EOM_constraints(t,x,aircraft,sim_options),[sim_options.t_start,sim_options.t_end],sim_options.ICS,sim_options.int_options);    
        else
            [t_out,X_out] = ode45(@(t,x)FLIGHT_SIM_EOM_constraints(t,x,aircraft,sim_options), sim_options.t, sim_options.ICS, sim_options.int_options);
        end
    end
    
end