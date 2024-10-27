% ground_hit.fcn terminates the simulation if the ground is hit.
%
% INPUTS:
%   t: time
%   x: state vector
%
% OUTPUTS:
%   value: 
%   isterminal:
%   direction:
%
% Sam Jaeger
% 2/12/2024
%   will eventually change to terminate_flight when I write my own ode45
%   solver.


function [value, isterminal, direction] = ground_hit(t,x)
    
    
    %%%%%% GROUND HIT %%%%%%%%%%%%%%%%%%%%%
    value = -x(9); % -z_f (altitude) is zero
    isterminal = 1; % yes, halt integration
    direction = 0; % cannot approach from negative numbers

    
end