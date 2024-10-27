%[ C_L, C_m_LE] = Coef_Vortex_Panel(vortex,coords,c,V_inf,alpha)
% Coef_Vortex_Panel.fcn computes the lift coefficient and moment
% coefficient from a given set vortex strengths and coordinates.
% Script uses Kutta Joukoski law to back out lift and pitching moment.
% Equations from Mechanics of flight 2nd edition, Pg. 36 (1.6.32, 1.6.33).
%
% Inputs:
%   vortex: vector of vortex panels
%   panel_length: vector of panel lengths
%   c: chord length
%   V_inf: free stream velocity
%   alpha: angle of attack (rad)
%
% Outputs:
%   C_L: lift coefficient
%   C_m_LE: moment coefficient at the Leading Edge
%
% Sam Jaeger
% Written: 1/12/2022
% Revised: 12/1/2022

function [C_L, C_m_LE] = Coef_Vortex_Panel(vortex,coords,c,V_inf,alpha)
    C_L = 0;
    C_m_LE = 0;
    
    N = length(vortex);
    
    l = zeros(N,1);
    for i=1:(N-1)
        dx = coords(i,1) - coords(i+1,1);
        dy = coords(i,2) - coords(i+1,2);
        l(i) = sqrt(dx^2 + dy^2);
    end
    x = coords(:,1);
    y = coords(:,2);
    
    
    for i=1:(N-1)
        C_L = C_L + l(i)*(vortex(i) + vortex(i+1));
        C_m_LE = C_m_LE + l(i)*( (2*x(i)*vortex(i)+x(i)*vortex(i+1)+x(i+1)*vortex(i)+2*x(i+1)*vortex(i+1))*cos(alpha) + (2*y(i)*vortex(i)+y(i)*vortex(i+1)+y(i+1)*vortex(i)+2*y(i+1)*vortex(i+1))*sin(alpha) );
    end
    
    C_L = C_L/(c*V_inf);
    C_m_LE = -C_m_LE/(3*V_inf*c^2);
end