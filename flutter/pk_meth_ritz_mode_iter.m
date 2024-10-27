% pk_meth_ritz_mode_iter.m script is an iterative p-k method using a Ritz 
% series of assumed shapes.
% This script will compute a frequency-damping-airspeed plots for a set of
% reduced frequencies. The calculations in this script are based on the
% equations of motion matrices in Section 5.6 of Hodges, Pierce
% Introduction to Structural Dynamics and Aeroelasticity. Additionally, the
% method for setting up and solving the eigenvalue problem was taken from
% Introduction to Unsteady Aerodynamics and Dynamic Aeroelasticity by
% Demasi in Chapter 36.5. This script differs from the pk_meth_ritz.m
% script in that it will compute all of the curves associated with the
% an inital guess mode natural frequency. The script will then attempt to
% plot only the damping and frequency of the target mode of interest.
% If modes are closely spaced, the algorithm often will pick the wrong mode
% which will appear as discontinuities in the V-omega and V-g plots.
% Additionally, this script is an update from pk_meth_ritz_k_mode_iter.m as
% the minimization function to pick the new k value is now a sqare instead 
% of an absolute value. Another update is that the inital guess for k uses
% the previously converged value of k for beyond the first freestream
% airspeed of interest. For the first freestream airspeed of interest, the
% natural frequency of the
%
% Sam Jaeger
% 8/28/2024

%% Inputs 

% number of bending modes
N_w = 3;

% number of torsion modes
N_theta = 3;

% Beam properties - in SI
props.m = 10;
props.l = 5;
props.e = -1/10;
props.a = -1/5;
props.b = 0.25;
props.r = sqrt(6/25);
props.EI = 1000;
props.GJ = 750;

% Goland Wing Beam Properties -in US
% props.m = 0.746;%10;
% props.l = 20;%5;
% props.e = -0.07;%-1/10;
% props.a = -0.17;%-1/5;
% props.b = 3;%0.25;
% props.r = 1.5;%sqrt(6/25);
% props.EI = 31.7e6*0.746;%1000;
% props.GJ = 1.23e6*1.943;%750;

% velocity space
u_start = 1;
u_end = 100; %*1.47; % to convert to ft/s to mi/hr

% number of flutter points to compute
%   in this case it is the number of u points
N_flutter_pts = 1000;

% altitude
alt = 0; % m

% Mach number guess
M_inf_guess = 0;

% reduced frequency convergence tolerance
tol = 0.001;

% Number of iterations to break the k convergence computation
% usually 10 or less
N_iter_break = 30;

%% Calculations

b = props.b;

% atmospheric properties
[rho_inf, T_inf, p_inf, a_inf, ~, ~, mu_inf] = ATMOS_1976(alt,'SI');

% mass, stiffness, coupling matrices
[M, K, A, omega_w, omega_theta] = assumed_shapes_mats(N_w, N_theta, props);

% damping matrix (if needed)
alpha = 0.00; % mass damping participation factor
beta = 0.000; % stiffness damping participation factor
C = alpha*M + beta*K; % Rayleigh Damping

% transform into modal space
[PHI,omega_n_2] = eig(K,M); % PHI is the modal matrix, omega_n_2 is a diagonal of natural frequencies squared
M_modal = PHI'*M*PHI; % should be diagonal
K_modal = PHI'*K*PHI; % should be omega_n_2
omega_n = [omega_w; omega_theta]; % vector of natural frequencies

% velocity 
u_inf = linspace(u_start, u_end, N_flutter_pts)';

% initalize
p = zeros((N_w+N_theta), length(u_inf), (N_w+N_theta));
g = zeros( length(u_inf), (N_w+N_theta));
g_full = zeros((N_w+N_theta), length(u_inf), (N_w+N_theta)); % all computed damping
omega2 = zeros( length(u_inf), (N_w+N_theta)); % real ones for plotting
omega = zeros((N_w+N_theta), length(u_inf), (N_w+N_theta)); % all computed frequencies
k_full = zeros(2*(N_w+N_theta), length(u_inf), (N_w+N_theta)); % all computed reduced frequencies
lambda = zeros(2*(N_w+N_theta), length(u_inf), (N_w+N_theta)); % eigenvalues
Mach = zeros( length(u_inf),1); 
k = zeros( length(u_inf),(N_w+N_theta)); % need to converge to a k value for each mode and velocity pair
iter = zeros( length(u_inf),(N_w+N_theta));
eps = zeros( length(u_inf),(N_w+N_theta)); %tolerance


for ii=1:length(u_inf) % loop over each value of u_inf
    for jj=1:(N_w + N_theta) % loop over each mode
        % inital guess for k...
        %   choice of natural frequency for inital guess affects results a lot
        if ii==1 % at the first airspeed of interest, guess the natural frequency of the target mode
            k(ii,jj) = (omega_n(jj))*b./u_inf(ii); 
        else % use the previous converged k as the inital guess
            k(ii,jj) = k(ii-1,jj);
        end

        eps(ii,jj) = 10; % error of k at each iteration, set something high to enter loop
        iter(ii,jj) = 1; % iteration counter
    
        while eps(ii,jj) > tol 
    
            [M_aero, C_aero, K_aero] = gen_aero_force_mats(k(ii,jj),N_w,N_theta,A,props);
            
            % assemble Aerodynamic influence coefficent matrix
            AIC = -(k(ii,jj)^2)*M_aero + (k(ii,jj)*b*1i)*C_aero + (b^2)*K_aero;
       
            % quadratic eigenvalue problem
            G = [zeros(N_w+N_theta), eye(N_w+N_theta); -M\(((b/u_inf(ii))^2)*K - rho_inf*AIC), -M\(b*C/u_inf(ii) )];
        
            lambda(:,ii,jj) = eig(G,'nobalance'); % compute eigenvalues
            % CAUTION!!! matlab will mix up the order of returned eigenvalues
            % could implement eigenshuffle algorithm, or just look at the plot
            %   https://mathworks.com/matlabcentral/fileexchange/22885-eigenshuffle
            % eigenshuffle.m probably won't work since each k needs to be
            % converged for each mode and airspeed combination.
    
            % check if guess is close 
            k_full(:,ii,jj) = imag(lambda(:,ii,jj)); % set of reduced frequencies
            
            % remove negative reduced frequency entries
            omega_counter = 1;
            for ll=1:length(k_full(:,ii,jj))
                if k_full(ll,ii,jj) >= 0
                    p(omega_counter,ii,jj) = lambda(ll,ii,jj);
                    omega(omega_counter,ii,jj) = k_full(ll,ii,jj)*u_inf(ii)/b;
    
                    g_full(omega_counter,ii,jj) = 2*real(lambda(ll,ii,jj))/k(ii,jj);
    
                    omega_counter = omega_counter + 1;
                end
            end
    
            % tolerance - minimize quadratic cost function
            [eps(ii,jj), min_index] = min( (omega(:,ii,jj)*b./u_inf(ii) - k(ii,jj)).^2 ); % find closest value to the guess reduced frequency
    
            % new reduced frequency
            k(ii,jj) = omega(min_index,ii,jj)*b./u_inf(ii); 
            
            % record omega and g values for the mode of interest
            omega2(ii,jj) = omega(min_index,ii,jj)^2;
            g(ii,jj) = g_full(min_index,ii,jj);

            iter(ii,jj) = iter(ii,jj)+1;
            % hopefully reduced freuqency converges in a normal amount of iterations - should be less than 10
            if iter(ii,jj) >= N_iter_break  
                disp(append(num2str(jj),' mode, k did not converge at U = ',num2str(u_inf(ii)) ))
                break
            end
        end
       
    end

    Mach(ii) = u_inf(ii)/a_inf;
end


%% Results and Plots

% find flutter point
%   check if sign of damping flips from negative to positive.
%   This isn't a good method since matlab will change the sorting of the
%   eigenvalues.
%   For now look at the plot when damping flips sign

for jj=1:(N_theta+N_w) % loop over i roots computed
        figure(369)
        subplot(2,1,1),plot(u_inf, sqrt(omega2(:,jj)),"."); hold on
        xlabel('$ U_\infty \left(\frac{m}{s}\right)$','FontSize',12,'Interpreter','latex')
        ylabel('$\omega  \left(\frac{rad}{s}\right)$','FontSize',12,'Interpreter','latex')
        %ylim([0, 20])
        xlim([0 max(u_inf)])
        grid on
        subplot(2,1,2),plot(u_inf, g(:,jj),"."); hold on
        xlabel('$ U_\infty \left(\frac{m}{s}\right)$','FontSize',12,'Interpreter','latex')
        ylabel('$g$','FontSize',12,'Interpreter','latex')
        grid on
        ylim([-.4, .2])
        xlim([0 max(u_inf)])
        hold on
    
        figure(370)
        subplot(2,1,1),plot(Mach, sqrt(omega2(:,jj)),"."); hold on
        xlabel('$ M_\infty $','FontSize',12,'Interpreter','latex')
        ylabel('$\omega  \left(\frac{rad}{s}\right)$','FontSize',12,'Interpreter','latex')
        %ylim([0, 20])
        xlim([0 0.8])
        grid on
        subplot(2,1,2),plot(Mach, g(:,jj),"."); hold on
        xlabel('$ M_\infty $','FontSize',12,'Interpreter','latex')
        ylabel('$g$','FontSize',12,'Interpreter','latex')
        grid on
        ylim([-.2, 0.2])
        xlim([0 0.8])
        hold on
        

    figure(400)
    plot(real(p(:,:,jj)),imag(p(:,:,jj)),'.')
    title('p-k method roots using Ritz series and Theodorsen Aerodynamics','FontSize',14,'Interpreter','latex')
    xlabel('Real p','FontSize',12,'Interpreter','latex')
    ylabel('Imag p','FontSize',12,'Interpreter','latex')
    grid on
    hold on

end
figure(369),subplot(2,1,1),title('p-k method using Ritz series and Theodorsen Aerodynamics','FontSize',14,'Interpreter','latex')
figure(370),subplot(2,1,1),title('p-k method using Ritz Series and Theodorsen Aerodynamics','FontSize',14,'Interpreter','latex')
hold off
