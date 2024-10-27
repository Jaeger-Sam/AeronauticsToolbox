% [coords] = NACA_4_dig_Vortex(m,p,TT,c,num_pts,cosine_cluster,disp_plot,close_open)  
%
% This function calculates the xy coordinates of a NACA 4 digit airfoil to
% be analyzed with a vortex panel script. The general outline of
% discritization can be found on Pg. 32 and 33 of Mechanics of Flight 2nd
% Ed. The script is set up so that the user can cluster points around the
% LE and TE using a cosine distribution. Additionally, a modified thickness
% can be inputted in order to close or leave the TE open. Generally better
% results are found by closing the TE and using a cosine clustering.
%
% Inputs:
%   m: maximum Camber in percent of chord
%   p: Location of camber in tenths
%   TT: Maximum thickness in % of the symmetric (uncambered) airfoil shape
%   c: Length of chord
%   num_pts: total number of points to define the airfoil (should be even)
%   cosine_cluster: ==1 to use a cosine cluster around the LE and TE
%   disp_plot: ==1 to display plot
%   close_open: ==0 for exact NACA 4 digit geometry (leaves TE open),
%               ==1 to close TE with a modified thickness. This gives a
%               better vorticity distribution around the TE.
%
% Outputs:
%   coords: Nx2 matrix of x and y locations of the airfoil geometry.
%           First column is x positions, second column is y position.
%           Starts on the upper surface from trailing edge and wraps around
%            to leading edge to the lower surface and then to the Trailing edge.
%
% Sam Jaeger
% Written: 1/7/2022
% Revised: 12/1/2022


function [coords] = NACA_4_dig_Vortex(m,p,TT,c,num_pts,cosine_cluster,disp_plot,close_open) 
    %% definitions
    if cosine_cluster == 1
        iseven  = rem(num_pts,2); 
        if iseven == 0 % even number of points, panel perpendicular to start of chord
            d_theta = pi / (real(num_pts/2) - 0.5);
            for i = 1:(num_pts/2)
                x(i) = c*0.5*(1 - cos(real(i)*d_theta - 0.5*d_theta));
            end
        elseif  iseven == 1 % odd number of points, coordinate on LE
            d_theta = pi / real(num_pts/2);
            x(1) = 0;
            for i = 2:floor(num_pts/2)
                x(i) = c*0.5*(1 - cos(real(i)*d_theta));
            end
        end
        x(end) = c;
    else
        x = linspace(0,c,floor(num_pts/2));  
    end
    
    camber_max = c*m/100; 
    t_max = c*TT/100; % maximum thickness
    x_c = x./c; %percent chord
    
    %% airfoil polynomial
    y_airfoil = zeros(length(x_c)-1,1);
    for i=1:(length(x_c)-1)
        %y_airfoil(i) = 5*t_max*(0.2969.*sqrt(x_c(i)) - 0.126.*(x_c(i)) - 0.3516.*(x_c(i)).^2 + 0.2843.*(x_c(i)).^3 - 0.1015.*(x_c(i)).^4);
        
        if close_open == 1
            % modified thickness to close trailing edge
            %y_airfoil(i)= tm*(2.980*m.sqrt(x) - 1.320*x - 3.286*x**2 + 2.441*x**3 - 0.815*x**4)
            y_airfoil(i) = 5*t_max*(2.980/10.*sqrt(x_c(i)) -  1.320/10.*(x_c(i)) - 3.286/10.*(x_c(i)).^2 + 2.441/10.*(x_c(i)).^3 - 0.815/10.*(x_c(i)).^4);
        else %exact coordinates
            y_airfoil(i) = 5*t_max*(0.2969.*sqrt(x_c(i)) - 0.126.*(x_c(i)) - 0.3516.*(x_c(i)).^2 + 0.2843.*(x_c(i)).^3 - 0.1015.*(x_c(i)).^4);
        end
    end
    y_airfoil = [y_airfoil; 0]; % make end point 0
    
    %% camber line
    y_c = zeros(length(x_c),1);
    for i = 1:length(x_c)   
        if x_c(i) < p/10
            y_c(i,1) = camber_max*((20/p)*x_c(i) - (100/(p^2))*(x_c(i)^2));
        else 
            y_c(i,1) = (c*m/(100 - 20*p + (p^2)))*(1 - 2*p/10 + (2*p/10)*x_c(i) - x_c(i)^2);
        end
    end
    
    %% airfoil surfaces
    
    upper = y_c + y_airfoil; % make upper surface of airfoil
    lower = y_c - y_airfoil; % make lower suface of airfoil
    x=x'; % make x position into a column vector
    

    %% Plot
    if disp_plot ==1
        figure,plot(x,y_airfoil,'-o',x,y_c,'-o')
        xlabel('chord position (x)')
        ylabel('thickness position (y)')
        legend('airfoil polynomial','chord polynomial')
        axis equal
        grid on

        figure,plot(x,upper,'o',x,lower,'o')
        xlabel('chord position (x)')
        ylabel('thickness position (y)')
        legend('upper surface of airfoil', 'lower surface of airfoil')
        axis equal
        grid on
    end
    
    %% coords 
    %lower surface then upper surface 
    %coords = [flip(x), flip(lower); x, upper];
    
    % starting at TE upper surface then lower surface
    coords = [ flip(x), flip(upper);  x, lower];
end