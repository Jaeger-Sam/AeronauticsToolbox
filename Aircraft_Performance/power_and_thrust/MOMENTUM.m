% MOMENTUM.fcn computs the relationships for momentum theoy of a propeller.
% This program is outlined on page 297 of Flight Mechanics by McCormick.
% All units must be in SI (meters, m/s, watts, newtons, etc.). Must have
% ATMOS.fcn in same working directory.
%
% INPUTS:
%   d_p: Propeller diameter (m)
%   V_inf: free stream velocity (m/s)
%   alt: altitude (m)
%
% VARIABLE INPUTS: 
%   Mutst be a vector with first entry being the value of interest (T_r or 
%   P_a) and the second being what the value is (thrust==1 or power==2)
%   
%   T_r: Thrust required in newtons [1] 
%   OR
%   P_a: Power available in watts [2]
%
% OUTPUTS:
%   w: induced velocity (meters / s)
%   eta_i: ideal efficiency (-)
%
% VARIABLE OUTPUTS:
%   P_i: Ideal power required in watts (if required thrust is specified)
%   OR
%   T_a: Thrust available in newtons (if power available is specified)
%
%
% Written by:
%   Sam Jaeger
%   1/31/2023

function [w,eta_i,Pi_Ta] = MOMENTUM(d_p,V_inf,alt,units,Tr_Pa)
    
    % error logic for power or thrust input
    if Tr_Pa(2) == 1 % thrust is specified
        T = Tr_Pa(1);
    elseif Tr_Pa(2) == 2 % power is specified
        P = Tr_Pa(1);
    else
        error('must specify thrust required == 1 or power available == 2')
    end

    A = pi*(d_p/2)^2;
    rho = ATMOS(alt,units);  

    % logic for power or thrust output
    if Tr_Pa(2) == 1 %thrust required is specifed, return power required 
       
        w = 0.5*(-V_inf + sqrt((V_inf^2)+ (2*T/rho/A) )); % equation 6.15, induced velocity
        eta_i = 1 / (1 + (w/V_inf)); %equation 6.19, efficiency
        Pi_Ta = T*(V_inf + w); % equation 6.15, power

    elseif Tr_Pa(2) == 2 % power available is specified, return thrust available
        % must iterate to get thrust

        % static values, inital guesses
        T0 = (P*sqrt(2*rho*A))^(2/3); % 6.17
        w0 = sqrt(T0/(2*rho*A)); % 6.16

        % use fzero, sub eq 6.12 into 6.15, solve for w
        w = fzero(@(w)(0.5*(-V_inf + sqrt((V_inf^2)+ (2/rho/A*(P/(V_inf+w) )) )) - w),w0);

        eta_i = 1 / (1 + (w/V_inf)); %equation 6.19, efficiency
        T = P/(V_inf + w); % equation 6.15, thrust
        Pi_Ta = T;
    end

end