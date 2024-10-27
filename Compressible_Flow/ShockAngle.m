% beta = ShockAngle(M1,delta,gamma)
%
% ShockAngle.m computes the shock turning angle beta given the upstream
% Mach number, the wedge angle, and ratio of specific heats (gamma) for a
% wedge. See Hammit, Murthy, "Approximate Solutions for Supersonic Flow 
% over Wedges and Cones" Journal of Aerospace Sciences, V 27, 1, 1960, pp. 
% 71-73. Returns only the weak shock solution.
%
% INPUTS:
%   M1: Upstream Mach number
%   delta: deflection angle (rad)
%   gamma: ratio of specific heats
%
% OUTPUTS:
%   beta: shock turning angle (rad)
%
% Sam Jaeger
% 2/19/2024

function beta = ShockAngle(M1,delta,gamma)
    % logic for checking inputs
    if gamma < 1
        error('gamma cannot be less than 1!')
    elseif M1 < 1
        error('Mach number must be supersonic!')
    elseif pi/2 < delta
        error('delta must less than 90 degrees!')
    elseif delta < 0
        error('Expansion wave! delta must be greater than zero!')
    end

    % define quantities for ease of use
    A = 1 + ((gamma-1)/2)*M1^2;
    C = 1 + ((gamma+1)/2)*M1^2;
    U = M1^2 - 1;
    R = (18*A*C*U + 27*(A^2) - (U*C)^2)/(4*U^3);
    
    % check if shock is detached
    delta_max = sqrt(acot(R/2 + sqrt((R^2 /4 ) + (A*C^3)/(U^3) ) ) );
    if delta > delta_max
        warning('DETATCHED SHOCK!')
        beta = NaN;
        return
    end

    % calculate beta
    D = cot(delta);
    B = -D*U;
    E = B^2 - 3*A*C;
    F = 2*E^(3/2);
    G = 9*A*B*C - 27*(A^2)*D - 2*B^3;

    phi = ( acos(G/F) + 4*pi)/3;
    if isreal(phi)
        beta = atan2(2*sqrt(E)*cos(phi)-B,  3*A);
    else
        error('DETACHED SHOCK!')
    end
end