% k method 2DOF script
%
% This script computes the flutter boundary using the k method described in
% 5.3.2 and 5.4.1 in Hodges, Pierce Intro. to Struc. Dynam. and
% Aeroelasticity. The physical parameters used for the 2DOF boundary are
% are described on page 186.
% See Page 193 and 194 for the steps in the k method loop
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

% Reduced frequency space
k_start = 0.001;
k_end = 3.0;

% Mach number guess
M_inf_guess = 0.0; 

%% k method loop

% define reduced frequency space
k = linspace(k_start, k_end, 1000)';


% initalize matrices
%   2 roots, length of k analysis points
Z = zeros(length(k),2);
omega_rat = zeros(length(k),2); % omega / omega_theta
g = zeros(length(k),2);
%clear k_F
for ii=1:length(k)
    
    % compute aero loads at each reduced freq.
    [l_h, l_theta, m_h, m_theta] = theodorsen(k(ii),a);

    % flutter det char. eqn.
    C2 = (sigma*mu*r)^2;
    C1 = -( (1+sigma^2)*(mu*r)^2 + mu*m_theta*sigma^2 + l_h*mu*r^2);
    C0 = mu*(r^2)*(mu + l_h) + m_theta*(mu + l_h) - (mu*x_theta + l_theta)*(mu*x_theta + m_h);

    % solve char. eqn.
    %   rows are a particular reduced frequency
    Z(ii,:) = roots([C2,C1,C0]);

    % find omega/omega_theta, g
    omega_rat(ii,:) = 1 ./ sqrt(real(Z(ii,:))); % omega / omega_theta
    %   only take positive solution
    g(ii,:) = imag(Z(ii,:)) ./ real(Z(ii,:)); 

    % This isnt working
%     % check if flutter point is hit
%     if exist("k_F","var")
%     elseif abs((g(ii,1))) < 0.01 || abs((g(ii,2))) < 0.01
%         k_F = k(ii);
%     end
end


%% Results and Plots

% find flutter point
%   check if sign of damping flips from negative to positive
for ii=1:(length(k)-1)
    if sign(g(ii,1)) ~= sign(g(ii+1,1))
        k_F = k(ii); % reduced frequency at flutter
        omega_F = omega_rat(ii,1); % frequency of flutter
        V_F = omega_rat(ii,1)/k_F; % reduced velocity at flutter
        break
    end
    
    if sign(g(ii,2)) ~= sign(g(ii+1,2))
        k_F = k(ii);
        omega_F = omega_rat(ii,2);
        V_F = omega_rat(ii,2)/k_F;
        break
    end
end

% Problem 5 in Chapter 5 says 
%   U_F = 2.170 b*omega_theta
%   omega_F = 0.6443 *omega_theta

% fake root locus
figure(70)
plot( real(Z), imag(Z), '.'); 
title('Root plot of $Z$','FontSize',14,'Interpreter','latex')
xlabel('$Re( Z )$','FontSize',12,'Interpreter','latex')
ylabel('$Im( Z )$','FontSize',12,'Interpreter','latex')
grid on; xlim([0 10]); ylim([-10 10])


% plot freq and damp vs U/b/omega_theta
for jj=1:2
    figure(60)
    subplot(2,1,1),plot(omega_rat(:,jj)./k, omega_rat(:,jj)); hold on
    xlabel('$ \frac{U}{b \omega_\theta}$','FontSize',12,'Interpreter','latex')
    ylabel('$\frac{\omega}{\omega_\theta}$','FontSize',12,'Interpreter','latex')
    ylim([0, 1.1])
    xlim([0 3])
    grid on
    subplot(2,1,2),plot(omega_rat(:,jj)./k, g(:,jj)); hold on
    xlabel('$ \frac{U}{b \omega_\theta}$','FontSize',12,'Interpreter','latex')
    ylabel('$g$','FontSize',12,'Interpreter','latex')
    grid on
    ylim([-0.9, 0.3])
    xlim([0 3])
    hold on
end
figure(60),subplot(2,1,1),title('Figures 5.18 and 5.19 in Hodges, Pierce Intro. to Struc. Dyn. and Aeroelas.','FontSize',14,'Interpreter','latex')
figure(60),subplot(2,1,1),xline(V_F,'-','$V_F$','Interpreter','latex')
figure(60),subplot(2,1,2),xline(V_F,'-','$V_F$','Interpreter','latex')
hold off

% these plots match Figures 5.18 and 5.19
