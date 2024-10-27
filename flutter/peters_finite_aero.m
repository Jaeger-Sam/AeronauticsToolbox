% peters_finite_aero.fcn computes the [A], {c}  matrices that define a
% finite set of aerodynamic states for a 2d airfoil. 
% See Secion 5.5.2 in Introduction to Structural Dynamics and 
% Aeroelasticity by Hodges, Pierce.
% See Peters, Karunamoorthy, Cao "Finite State Induced Flow Models Part I:
% Two-Dimensional Thin Airfoil", in Journal of Aircraft 1995
%
% INPUTS:
%   N_aero_states: number of aerodynamic states. Having >10 will make the 
%       problem ill conditioned 
%
% OUTPUTS:
%   A: N x N matrix defining the aero states dynamics
%   c: N x 1 vector defining the aero states forcing from airfoil motion
%
% Sam Jaeger
% 7/31/2024

function [A, c] = peters_finite_aero(N_aero_states)
    N = N_aero_states;
    
    D = zeros(N);
    d = zeros(N,1);
    b = zeros(N,1);
    c = zeros(N,1);

    for n = 1:N
        % populate D
        for m = 1:N
            if n == m+1
                D(n,m) = 1/2/n;
            elseif n == m-1
                D(n,m) = -1/2/n;
            end
        end

        % populate b
        if n ~= N
            b(n) = ((-1)^(n-1)) *factorial(N+n-1)/factorial(N-n-1)/(factorial(n)^2);
        elseif n == N
            b(n) = ((-1)^(n-1));
        end

        % populate d
        if n == 1
            d(n) = 0.5;
        end
        
        % populate c
        c(n) = 2/n;

    end

    A = D + d*b' + c*d' + 0.5*c*b';
end