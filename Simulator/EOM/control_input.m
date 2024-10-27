% This function retuns a vector of control inputs given the state vector,
% nominal y output, and C matrix where {y} = [C]*{x}. ([D] is assumed to be
% zero). Also gain matrix K needs to be inputed.
%
% INPUTS:
%   x: state vector
%   yb: Nonminal y output 
%   C: matrix mapping state vector to sensor output
%   K: gain matrix
%   max_deflect: matrix of saturation states (rad)
%
% OUTPUTS:
%   u: vector of control inputs
%
% Sam Jaeger
% 1/3/2024
%   Revised: 1/23/2024
%   Revised: 3/15/2024

function u = control_input(x,yb,C,K,max_deflect,t)
    y = C*x;
    %u = K*(yb([1,3,5]) - y);
    %u = K*y;

    % for LQR controller developed...
    %yb = [7.7624e+03;1.9354e+03;0;8000]; % hard code ybar in for now.
    %u = -K*([x(1);x(3);x(5);norm(x(1:3))] - yb);
    u=0;
    
    u = [0;u;0;0];

    % Deflect 5 deg for the first 20 seconds
%     if t < 20
%         u(2) = (-4+5)*pi/180;
%     else
%         u(2)= -4*pi/180;
%     end

    % Saturation Limits
    for ii=1:length(u)
        if u(ii) > max_deflect(ii,1)
            u(ii) = max_deflect(ii,1);
        elseif u(ii) < max_deflect(ii,2)
            u(ii) = max_deflect(ii,2);
        end
    end
end