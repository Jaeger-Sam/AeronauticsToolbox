clear all
% close all

% Example numbers from:
% https://ftp.unpad.ac.id/orari/library/library-non-ict/aero/docs/LIFTING%20LINE%20THEORY%20Tutorial%20Example.pdf
% 8 linearly distrbuted nodes
% C_L =0.2325
% C_D_i= 0.002253
% gave e=0.9546



% Lifting Line Checkout
N = 100;

C_L_alpha_2d = 6;
alpha_root = -8*(pi/180);
alpha_L0_root = -2*(pi/180);
c_root = 1.875;
c_tip = 3.125;
b=20;
V_inf = 100;
washout = 2.9*pi/180; %2.9

% cosine distribution
theta=linspace(0,pi,N);
for i=1:(N/2)
    c(i,1)=(c_root-c_tip)*cos(theta(i)) + c_tip;
end
c=[c;flip(c)]';
%c=[linspace(c_root,c_tip,N/2),flip(linspace(c_root,c_tip,N/2))]';
%c = 1*ones(N,1);
 
 


%A = Lifting_Line_Coef(alpha,C_L_alpha_2d,alpha_L0,b,c,N);

[C_L,C_D_i,oswald_efficiency_factor,z,gamma] = Lifting_Line(alpha_root,V_inf,washout,alpha_L0_root,C_L_alpha_2d,c,b,N);

C_L
C_D_i
oswald_efficiency_factor  

figure,plot(z,gamma)
xlabel('span position')
ylabel('vorticity')