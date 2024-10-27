% code to compute equilibrium properties of air given pressure and enthalpy
% effects of ionization and electronic excitation are not included.
%
% INPUTS:
%   p1: pressure in Pa
%   h1stag: enthalpy in MJ/kg
%   disp_results: logical to print results or not
%
% OUTPUTS:
%   Z: 
%   geff:
%   hchecm: enthalpy J/kg
%   cs: N2,O2,NO,N,O mass fractions = 
%   t1: temperature K
%   r1: density, kg/m^3
%
% Adopted from Graham Candler AEM 5247 class code
% Sam Jaeger 
% 2/4/2024

function [Z,geff,hchem,cs,t1,r1] = mollier(p1,h1stag,disp_results)
      %clear

      r = 8314.3;
      b0m = 0.00;
      
      %p1 = input('pressure in Pa ');
      %h1stag = input('enthalpy in MJ/kg ');
      h1stag = h1stag*1.0e+06;

      a(1,1) = 0.97919118E+01;
      a(1,2) =-0.11280844E+02;
      a(1,3) =-0.79705999E-02;
      a(1,4) = 0.29358337E-03;
      a(1,5) =-0.31783590E-05;
      a(2,1) = 0.96749680E+01;
      a(2,2) =-0.58273255E+01;
      a(2,3) =-0.13627750E-01;
      a(2,4) = 0.45994475E-03;
      a(2,5) =-0.47994368E-05;
      a(3,1) = 0.82209491E+01;
      a(3,2) =-0.74585635E+01;
      a(3,3) =-0.10356698E-01;
      a(3,4) = 0.35911437E-03;
      a(3,5) =-0.37928553E-05;

      m(1)=28.016;
      m(2)=32.0;
      m(3)=30.008;
      m(4)=14.008;
      m(5)=16.0;

      thv(1) = 3395.0;
      thv(2) = 2239.0;
      thv(3) = 2817.0;
       
      h0(1) = 0.0;
      h0(2) = 0.0;
      h0(3) = 2.996120d6;
      h0(4) = 3.362160d7;
      h0(5) = 1.542000d7;

      cN2 = 0.767083;
      cO2 = 1.0-cN2;
      
      gcon = (cN2/m(1)+cO2/m(2))*r;

      beta = (cN2/m(1))/(cO2/m(2));

      t1 = h1stag/1000.0;
      h = 0.0;

      while abs(h-h1stag) > 0.1

      r1 = p1/(t1*r*(cN2/m(1)+cO2/m(2))+p1*b0m);
      xs(1) = r1*cN2/m(1);
      xs(2) = r1*cO2/m(2);
      xs(3) = 0.0;
      xs(4) = 0.0;
      xs(5) = 0.0;

      pp = 0.0;
      
      while abs(pp-p1) > 0.001

          z = 10000.0/t1;
          for n=1:3
            rk(n)= exp(a(n,1)+a(n,2)*z+a(n,3)*z*z+a(n,4)*z*z*z+a(n,5)*z*z*z*z);
          end
          
          for it=1:320
            f1 = 0.0;
            for n=1:5
                f1 = f1 + xs(n)*m(n);
            end
            f1 = f1 - r1;
    
            f2 = 2.0*xs(1)-2.0*beta*xs(2)+(1.0-beta)*xs(3)+xs(4)-beta*xs(5);
    
            df1d1 = m(1) + 0.5*(m(3)*xs(3)+m(4)*xs(4))/xs(1); 
            df1d2 = m(2) + 0.5*(m(3)*xs(3)+m(5)*xs(5))/xs(2); 
            df2d1 = 2.0 + 0.5*(xs(4)+(1.0-beta)*xs(3))/xs(1);
            df2d2 =-2.0*beta - 0.5*(beta*xs(5)+(1.0-beta)*xs(3))/xs(2);
    
            det = df1d1*df2d2-df1d2*df2d1;
            xs(1) = xs(1) - 0.5*( df2d2*f1-df1d2*f2)/det;
            xs(2) = xs(2) - 0.5*(-df2d1*f1+df1d1*f2)/det;
            xs(4) = sqrt(rk(1)*xs(1));
            xs(5) = sqrt(rk(2)*xs(2));
            xs(3) = xs(4)*xs(5)/rk(3);
          end
    
          pp = 0.0;
          r1 = 0.0;
          for n=1:5
            pp = pp + r*t1*xs(n);
            r1 = r1 + xs(n)*m(n);
          end
    
          pp = pp/(1.0-r1*b0m);
    
          r1 = r1*p1/pp;

      end
      
      for n=1:5
        cs(n) = xs(n)*m(n)/r1;
      end

      h = 0.0;
      for n=1:3
          h = h + cs(n)*(2.5*r*t1/m(n) + r*t1/m(n)/(1.0-r1*b0m) + h0(n)...
              +r/m(n)*thv(n)/(exp(thv(n)/t1)-1.0));
      end
      for n=4:5
        h = h + cs(n)*(1.5*r*t1/m(n) + r*t1/m(n)/(1.0-r1*b0m) + h0(n));
      end
 
      cp = h/t1;
      t1 = t1 - 0.25*(h-h1stag)/cp;

      end

     
      
      geff = (r1*h1stag/p1)/((r1*h1stag/p1)-1.0);
      Z = p1/(r1*gcon*t1);
      
      hchem = cs(1)*h0(1)+cs(2)*h0(2)+cs(3)*h0(3)+cs(4)*h0(4)+cs(5)*h0(5);

      if disp_results == true
          disp(sprintf('FINAL RESULTS FOR CONDITIONS:'));
          disp(sprintf(' h  = %10.2f J/kg',h1stag));
          disp(sprintf(' p  = %10.2f Pa',p1));
          disp(sprintf(' T  = %10.2f K',t1));
          disp(sprintf('rho = %10.5f kg/m^3',r1));
    
          disp(sprintf(' Z  = %10.5f',Z));
          disp(sprintf('geff= %10.5f',geff));
    
          disp(sprintf('hchem= %10.2f J/kg',hchem));
    
          disp(sprintf(' '));
          disp(sprintf('N2,O2,NO,N,O mass fractions = '));
          disp(sprintf('%9.5f',cs(1)));
          disp(sprintf('%9.5f',cs(2)));
          disp(sprintf('%9.5f',cs(3)));
          disp(sprintf('%9.5f',cs(4)));
          disp(sprintf('%9.5f',cs(5)));
      end
  end