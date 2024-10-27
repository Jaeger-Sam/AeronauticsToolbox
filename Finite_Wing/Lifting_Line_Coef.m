% Lifting_Line_Coef.fcn calculates the fouier coefficients for a wing with
% no sweep, no dihedral, and no sideslip.
% The wing can have varying cross section, which must be denoted by the
% alpha_L0 and C_L_alpha_2d vectors. Additionally the wing can have washout (twist) which
% must be denoted by alpha vector. Finally, the wing can have varying chord
% length which is denoted by c. All of these quanitites must be the same
% length as the number of horseshoe vorticies used.
%
% Inputs:
%   alpha: Nx1 vector of aerodynamic angle of attacks with respect to the span position 
%   C_L_alpha_2d: Nx1 vector of derivative of lift coefficient with respect to angle of attack for each section of wing.
%   alpha_L0: Nx1 vector of zero lift angle of attack with respect to the span position
%   b: scalar span of wing 
%   c: Nx1 vector of chord lengs with respect to the span position
%   theta: Nx1 vector of wpanwise locations in the theta coordinate system
%   num_point: (==N) integer scalar. number of horseshoe vorticies to model the lifting line. Should be an even number
%
% Outputs:
%   A: (N-2)x1 vector of fouier coefficients used to compute vortex (Gamma)
%       distribution on the lifting line. This vector is two less than the
%       total number of horseshoe segments becuase of the two wingtips. Even
%       indices are zero.
%
% Sam Jaeger
% 1/18/2021



function [A] = Lifting_Line_Coef(alpha,C_L_alpha_2d,alpha_L0,b,c,theta,num_point)
    N = num_point;

    %% error logic
    alpha_length = length(alpha);
    C_L_alpha_length = length(C_L_alpha_2d);
    alpha_L0_length = length(alpha_L0);
    c_length = length(c);
    
    if alpha_length ~= N || alpha_L0_length ~= N || c_length ~= N || C_L_alpha_length ~= N
        error('alpha, C_L_alpha, alpha_L0, c must be the Nx1 vectors')
    end


    %% Fill in M matrix of fourier summation
    M=zeros(N-2,N-2);
    % do not include first or last indicies (wingtips). We know that
    % vorticity (gamma) must equal zero at the wingtips. 
    for n=2:(N-1) %series values (columns)
        for i=2:(N-1) %theta locations (rows)
             M(i-1,n-1) = sin((n-1)*theta(i))*( (4*b/C_L_alpha_2d(i)/c(i)) + ((n-1)/sin(theta(i))) );
        end
    end
    
    %% Fill in dalapha vector
    dalpha=zeros(N-2,1);
    %do not include wingtips for reasons above
    for i=1:(N-2)
        dalpha(i) = alpha(i+1) - alpha_L0(i+1);
    end
    

    %% calculate A
    A = M\dalpha;


end