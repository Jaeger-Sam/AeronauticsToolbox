% k_meth_ritz.m script is a k method using a Ritz series of assumed shapes.
% This script will compute a frequency-damping-airspeed plots for a set of
% reduced frequencies. The calculations in this script are based on the
% equations of motion matrices in Section 5.6 of Hodges, Pierce
% Introduction to Structural Dynamics and Aeroelasticity. Additionally, the
% method for setting up and solving the eigenvalue problem was taken from
% Introduction to Unsteady Aerodynamics and Dynamic Aeroelasticity by
% Demasi in Chapter 36.3. Note, that because the k method you have to
% divide by k, you cannot solve at 0 reduced frequency. Thus, divergence
% cannot be computed.
%
% Sam Jaeger
% 8/1/2024

%% Inputs 

% number of bending modes
N_w = 2;

% number of torsion modes
N_theta = 2;

% Beam properties
props.m = 10;
props.l = 5;
props.e = -1/10;
props.a = -1/5;
props.b = 0.25;
props.r = sqrt(6/25);
props.EI = 1000;
props.GJ = 750;

% Goland Wind beam properties - in US
% props.m = 0.746;%10;
% props.l = 20;%5;
% props.e = -0.07;%-1/10;
% props.a = -0.17;%-1/5;
% props.b = 3;%0.25;
% props.r = 1.5;%sqrt(6/25);
% props.EI = 31.7e6*0.746;%1000;
% props.GJ = 1.23e6*1.943;%750;
% props.rho_Ip = 20;

% Reduced frequency space
k_start = 0.001; % cannot start at 0 using the k method
k_end = 8.0;

% number of flutter points to compute
%   in this case it is the number of k points
N_flutter_pts = 1000;

% altitude
alt = 0; % m

% Mach number guess
M_inf_guess = 0;



%% Calculations

b = props.b;

% atmospheric properties
[rho_inf, T_inf, p_inf, a_inf, ~, ~, mu_inf] = ATMOS_1976(alt,'SI');

% mass, stiffness, coupling matrices
[M, K, A, omega_w, omega_theta] = assumed_shapes_mats(N_w, N_theta, props);

% damping matrix (if needed)
%C = zeros(N_w+N_theta); % for now zero damping
alpha = 0.00; % mass damping participation factor
beta = 0.000; % stiffness damping participation factor
C = alpha*M + beta*K; % Rayleigh Damping

% transform into modal space
[PHI,omega_n_2] = eig(K,M); % PHI is the modal matrix, omega_n_2 is a diagonal of natural frequencies squared
M_modal = PHI'*M*PHI; % should be diagonal
K_modal = PHI'*K*PHI; % should be omega_n_2

% reduced frequency 
k = linspace(k_start, k_end, N_flutter_pts)';

% initalize
g = zeros(2*(N_w+N_theta), length(k));
omega2 = zeros(2*(N_w+N_theta), length(k));
U_inf = zeros(2*(N_w+N_theta), length(k));
Mach = zeros(2*(N_w+N_theta), length(k));

for ii=1:length(k)
    [M_aero, C_aero, K_aero] = gen_aero_force_mats(k(ii),N_w,N_theta,A,props);
    
    % assemble Aerodynamic influence coefficent matrix
    AIC = -2*((k(ii)/b)^2)*M_aero + 2*(k(ii)/b)*1i*C_aero + 2*K_aero;

    % new aeroelastic system mass matrix for quadratic eigenvalue problem
    D = M + 0.5*rho_inf*((b/k(ii))^2)*AIC; % 36.42

    % quadratic eigenvalue problem
    B = [zeros(N_w+N_theta), eye(N_w+N_theta); -D\K, -D\C]; %36.53

    lambda = eig(B,'nobalance'); % compute eigenvalues
    %   nobalance seems to help but not eliminate the sorting problem
    % CAUTION!!! matlab will mix up the order of returned eigenvalues
    % could implement eigenshuffle algorithm, or just look at the plot
    %   https://mathworks.com/matlabcentral/fileexchange/22885-eigenshuffle
    
    % two eigenvalues for each mode
    % damping and frequency should be the same for each mode if there is
    % zero structural damping.
    g(:,ii) = - imag(lambda.^2)./real(lambda.^2); % damping
    omega2(:,ii) = - (real(lambda.^2).^2 + imag(lambda.^2).^2)./real(lambda.^2); % frequencies squared
    U_inf(:,ii) = b*sqrt(omega2(:,ii))./k(ii); % freestream at each k point
    Mach(:,ii) = U_inf(:,ii)/a_inf;
end


%% Results and Plots

% find flutter point
%   check if sign of damping flips from negative to positive.
%   This isn't a good method since matlab will change the sorting of the
%   eigenvalues.
%   For now look at the plot when damping flips sign

for jj=1:2*(N_theta+N_w)
    figure(69)
    subplot(2,1,1),plot(U_inf(jj,:), sqrt(omega2(jj,:)),"."); hold on
    xlabel('$ U_\infty \left(\frac{m}{s}\right)$','FontSize',12,'Interpreter','latex')
    ylabel('$\omega  \left(\frac{rad}{s}\right)$','FontSize',12,'Interpreter','latex')
    %ylim([0, 100])
    xlim([0 100])
    grid on
    subplot(2,1,2),plot(U_inf(jj,:), g(jj,:),"."); hold on
    xlabel('$ U_\infty \left(\frac{m}{s}\right)$','FontSize',12,'Interpreter','latex')
    ylabel('$g$','FontSize',12,'Interpreter','latex')
    grid on
    ylim([-.2, .2])
    xlim([0 100])
    hold on

    figure(70)
    subplot(2,1,1),plot(Mach(jj,:), sqrt(omega2(jj,:)),"."); hold on
    xlabel('$ M_\infty $','FontSize',12,'Interpreter','latex')
    ylabel('$\omega  \left(\frac{rad}{s}\right)$','FontSize',12,'Interpreter','latex')
    %ylim([0, 20])
    %xlim([0 0.8])
    grid on
    subplot(2,1,2),plot(Mach(jj,:), g(jj,:),"."); hold on
    xlabel('$ M_\infty $','FontSize',12,'Interpreter','latex')
    ylabel('$g$','FontSize',12,'Interpreter','latex')
    grid on
    ylim([-.2, 0.2])
    xlim([0 0.8])
    hold on
    
end
figure(69),subplot(2,1,1),title('k method using Ritz series and Theodorsen Aerodynamics','FontSize',14,'Interpreter','latex')
figure(70),subplot(2,1,1),title('k method using Ritz Series and Theodorsen Aerodynamics','FontSize',14,'Interpreter','latex')
hold off