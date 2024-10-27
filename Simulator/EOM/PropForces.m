% PropForces.fcn computes the propulsion forces and moments given the state
% of the aircraft, control_input, and propulsion propertes.
%
% INPUTS:
%   x: state vector of aircraft
%   control_vec: vector of control inputs
%   propulsion: data structure of propulsion properties: 
%       aircraft.propulsion from initalize_sim.fcn
%
% OUTPUTS:
%   F_P: vector of propulsive forces in body fixed x,y,z coords
%   M_P: vector of moments in body fixed x,y,z coords
%
% Written By:
% Sam Jaeger
% 10/20/2023
%   Revised: 1/5/2024
%   Revised: 1/15/2024, added state vector input

function [F_P, M_P] = PropForces(x, control_vec, propulsion)
    delta_T = control_vec(1);

    % for now jet aircraft:
    F_P(1) = propulsion.max_thrust*delta_T;
    F_P(2) = 0;
    F_P(3) = 0;

    M_P = [0;0;0]; % no moments applied from thrust

end

