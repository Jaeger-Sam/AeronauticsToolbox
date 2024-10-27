% Lifting_Line.fcn calculates lift and induced drag coefficients for a
% finite wing. Additionally it returns the oswald efficiency factor (e),
% spanwise locations (z) and lifting line circulation distribution (gamma).
% No sideslip, dihedral, or sweep is allowed. Airfoil section is assumed to
% be constant throughout the span which is defined by C_L_alpha_2d and
% alpha_L0_root quanitites. This function calls Lifting_Line_Coef.fcn
%
% Inputs:
%   alpha_root: angle of attack of the root
%   V_inf: freestream velocity
%   washout: (rad) amount of twist in the wing (assumes linear distribution). (+) value will make tip at lower angle than the root.
%   alpha_L0_root: the zero lift angle of attack at the root of the wing
%   C_L_alpha_2d: scalar derivative of lift coefficient with respect to angle of attack for each section of wing.
%   c: Nx1 vector of chord lengths with respect to the span position
%   b: scalar span of wing
%   num_point: (==N) integer scalar. number of horseshoe vorticies to model the lifting line. This should be an even number so that this is no
%               control point on the center of the wing.
%   plot_fourier: plots the fouier coefficients normalized by angle of attack
%
% Outputs:
%   C_L: lift coefficient for entire wing
%   C_D_i: induced drag coefficient for the entire wing
%   oswald_efficiency_factor: e for induced drag calculations
%   z: Nx1 vector of spanwise locations where each horshoe vortex is computed. This is cosine clusted (more nodes by the wingtips). + is left wing - is right wing
%   gamma: Nx1 vector of circulation along the lifting line
%
% Sam Jaeger
% Revised: 10/30/2022
% Revised: 12/20/2022

function [C_L,C_D_i,oswald_efficiency_factor,z,gamma] = Lifting_Line(alpha_root,V_inf,washout,alpha_L0_root,C_L_alpha_2d,c,b,num_point)
    %% Calculation of fundamental quanitites
    N = num_point;

    %% Divide up spanwise locations - cosine cluster
    theta = linspace(0,pi,N);
    z=zeros(N,1);
    for i=1:N
        z(i,1) = (b/2)*cos(theta(i));
    end

    alpha = alpha_root - washout*abs(cos(theta));
    alpha_L0 = alpha_L0_root*ones(N,1);
    C_L_alpha_2d = C_L_alpha_2d*ones(N,1);

    %% Calc Fourier Coef
    A = Lifting_Line_Coef(alpha,C_L_alpha_2d,alpha_L0,b,c,theta,N);


    
    %% Calc vorticity distribution
    gamma=zeros(N,1);
    for i=2:(N-1)
        for n=2:(N-1)
            gamma(i) = 2*b*V_inf*A(n-1)*sin((n-1)*theta(i)) + gamma(i);
        end
    end

    %% Calculate Wing area
    d_area = zeros(N,1);
    dz=zeros(N,1);
    for i=1:N-1 % area of ith section
        dz(i)=z(i+1)-z(i);
        d_area(i) = abs((dz(i))*c(i));
    end
    S = sum(d_area); %area of wing
    AR = b^2/S ;% aspect Ratio

    %% Calculate lift and drag coefficients of wing

    C_L = A(1)*pi*AR;

    C_L_gamma_summation = 2*sum(gamma.*abs(dz))/V_inf/S;

    C_d_i_sum = 0;
    for n=2:(N-2)
        C_d_i_sum = C_d_i_sum + n*(A(n)/A(1))^2;
    end

    C_D_i = (pi*AR*A(1)^2)*(1 + C_d_i_sum);

    oswald_efficiency_factor = (C_L^2)/(pi*AR*C_D_i);

    %%
    format long
    AoA = alpha_root*180/pi
    theta'
    a = A/alpha_root

end