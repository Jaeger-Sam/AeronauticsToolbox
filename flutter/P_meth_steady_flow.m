% 2 DOF, steady aerodynamics p method script.
%
% This script will solve for the p roots of a 2DOF flutter system with
% using steay flow theory (C_L = 2*pi*alpha, M_(c/4) = 0). Details are
% given in Section 5.2 of Introduction to Structural Dynamics and 
% Aeroelasticity by Hodges and Pierce. Specifically see pages 185 & 186.
%
% Sam Jaeger
% 7/30/2024

%% Define inputs
% see page 186 for values

%  nondimensional percent distance from the half chord locating the shear
%  center (point "P"). 
%  When  a = -1, the shear center is at the LE. 
%  When  a = +1, the shear center is at the TE.
%  When  a =  0, the shear center is at the half chord (c/2 = b)
%       -1 <= a <= 1
a = -1/5;

%  nondimensional percent distance from the half chord locating the center
%  of mass (point "C"). 
%  When  e = -1, the mass center is at the LE. 
%  When  e = +1, the mass center is at the TE.
%  When  e =  0, the mass center is at the half chord (c/2 = b)
%       -1 <= e <= 1
e = -1/10;

% static unbalance
%   Nondimensional distance between center of mass and shear center.
%   Positive when center of mass is further aft than the shear center.
%   Negative when the shear center is further aft than the center of mass.
x_theta = e - a;

% mass ratio
%   mu = m / (rho_inf * pi * b^2)
mu = 20;

% radius of gyration about shear center
%   r^2 = I_p / (m * b^2)
%   I_p = I_C + m*(b^2)*(x_theta^2)
%       I_C is the moment of inertia about the center of mass
r = sqrt(6/25);

% ratio of natural frequencies of the springs attached at point P
%   sigma = omega_h / omega_theta
%       omega_h = sqrt( k_h / m )
%               = heave or vertical spring nat freq
%       omega_theta = sqrt( k_theta / I_p )
%                   = pitch or rotational spring nat freq
sigma = 2/5;

% reduced velocity analysis start and end values
%   V = U / (b omega_theta)
V_red_start = 0.0001; % cannot start at 0, since there are 1/V terms!
V_red_end = 3;

%% Assemble Polynomial & Solve

% define reduce velocity space
V = linspace(V_red_start, V_red_end, 300);

% initalize solved variables
p = zeros(length(V),4); % should have 4 roots for binary system
Gamma_omega_t = zeros(length(V),4);
Omega_omega_t = zeros(length(V),4);
for ii=1:length(V)
    C4 = r^2 - x_theta^2;
    C3 = 0;
    C2 = (sigma*r/V(ii))^2 + (r/V(ii))^2 - (2/mu)*(a + 0.5) - x_theta*2/mu;
    C1 = 0;
    C0 = (sigma*r/V(ii)^2)^2 - (2/mu)*(a + 0.5)*(sigma/V(ii))^2;
    
    % rows are particular reduced velocity
    % columns are particular root
    % there should be 2 complex conjuate pairs of roots
    p(ii,:) = roots([C4,C3,C2,C1,C0]);

    % re and im of evals and multiply them by reduced velocity to get
    % Gamma/omega_theta and Omega/omega_theta
    %   See eqn 5.34
    Gamma_omega_t(ii,:) = real(p(ii,:))*V(ii);
    % only really care about positive values since these are complex conjugate pairs
    Omega_omega_t(ii,:) = imag(p(ii,:))*V(ii); 
end



%% Print Solution and Plots

% divergence speed (reduced velocity)
%   when p = 0, solve for V ==> C0 = 0
V_D = r*sqrt(mu / (1 + 2*a))

% book says V_D = 2.828 for the given values

% find flutter reduced velocity
clear V_F Omega_F_omega_t
for ii=1:length(V)
    for jj=1:4
        if Gamma_omega_t(ii,jj) > 0.0001
            V_F = V(ii)
            Omega_F_omega_t = mean(abs(Omega_omega_t(ii,1))) % positive values should converge to the same number
            break
        end
    end
    if exist("V_F","var")
        break
    end
end

% books says V_F = 1.843, Omega_F / omega_theta = 0.5568

% plot freq and damp plots
for jj=1:4
    figure(53)
    subplot(2,1,1),plot(V, Omega_omega_t(:,jj)); hold on
    xlabel('$\frac{U}{b \omega_\theta}$','FontSize',12,'Interpreter','latex')
    ylabel('$\frac{\Omega}{\omega_\theta}$','FontSize',12,'Interpreter','latex')
    ylim([0, 1.1])
    grid on
    subplot(2,1,2),plot(V, Gamma_omega_t(:,jj)); hold on
    xlabel('$\frac{U}{b \omega_\theta}$','FontSize',12,'Interpreter','latex')
    ylabel('$\frac{\Gamma}{\omega_\theta}$','FontSize',12,'Interpreter','latex')
    grid on
    ylim([-0.6, 0.6])
    hold on
end
figure(53),subplot(2,1,1),title('Figures 5.3 and 5.4 in Hodges, Pierce Intro. to Struc. Dyn. and Aeroelas.','FontSize',14,'Interpreter','latex')
figure(53),subplot(2,1,1),xline(V_F,'-','$V_F$','Interpreter','latex')
figure(53),subplot(2,1,2),xline(V_F,'-','$V_F$','Interpreter','latex')
hold off

% all of these numbers and plots match the book