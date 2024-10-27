% [M2,P0_rat,T_rat,P_rat,rho_rat] = normal_shock(M1,gamma)
%
% normal_shock.fcn computes the properties across a normal shock given an 
% upstreem Mach and ratio of specific heats. Assumes a perfect, non
% reacting gas.
%
% INPUTS:
%   M1: Upstream Mach number
%   gamma: ratio of specific heats c_p/c_v
%
% OUTPUTS:
%   M2: Mach after shock
%   P0_rat: P02/P01
%   T_rat: T2/T1
%   P_rat: P2/P1
%   rho_rat: rho2/rho1
%
% Sam Jaeger
% 2/8/2024

function [M2,P0_rat,T_rat,P_rat,rho_rat] = normal_shock(M1,gamma)
    
    if M1 < 1
        error('Upstream Mach number must be supersonic!')
    elseif gamma < 1
        error('c_p/c_v must be greater than 1!')
    end

    g= gamma; % to simplify code
    % t1, t2, etc. are just to break up terms in normal shock expression
    % for easy coding

    t1 = 1 + 0.5.*(g-1).*M1.^2;
    t2 = (g.*M1.^2)-(g-1)*0.5;
    M2 = sqrt(t1./t2);

    t1=1+ 0.5*(g-1).*M1.^2;
    t2=1+ 0.5*(g-1).*M2.^2;
    T_rat = t1./t2;

    t1 = 1+ g.*M1.^2;
    t2 = 1+ g.*M2.^2;
    P_rat = t1./t2;

    rho_rat = P_rat./T_rat; % assuming perfect gas

    t1=(g+1).*(M1.^2);
    t2=2+ ((g-1).*(M1.^2));
    t3=g./(g-1);
    t4=(2.*g.*(M1.^2)) - (g-1);
    t5=1./(g-1);
    P0_rat = ((t1./t2).^t3).*((g+1)./t4).^t5;
end