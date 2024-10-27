clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NACA 4 digit inputs
m = 5;
p = 4;
TT = 06;

% chord length
c = 1;

% panel options for function generator
num_pts = 100;
cosine_cluster = 1;
disp_plot = 1;
close_open = 1;

% alpha sweep options
alpha_deg = linspace(-5,5,11);

% ratio of specific heats
gamma = 1.4;

%compressibility correction rule
comp_corr_rule = 2; % ==1 for Prandtl-Gaulert, ==2 for Karman-Tisen, ==3 for Laitone

% Panel options
panel_opts.V_inf=100;
panel_opts.Mach=0.4;
panel_opts.comp_corr_rule = comp_corr_rule;
panel_opts.gamma = gamma;
panel_opts.offset=0.00001; 
panel_opts.kj=true;
panel_opts.cp=true;
panel_opts.critical_mach = true;
panel_opts.sweep_updates=true;
panel_opts.alpha_plots=true;
panel_opts.fit_data = true;
panel_opts.alpha_fit_start = -2;
panel_opts.alpha_fit_end = 2;
panel_opts.cp_plots = true;


%% Function Calls

% get geometry
coords = NACA_4_dig_Vortex(m,p,TT,c,num_pts,cosine_cluster,disp_plot,close_open) ;

% send to panel code
data_out = panel_alpha_sweep(alpha_deg, coords, panel_opts);