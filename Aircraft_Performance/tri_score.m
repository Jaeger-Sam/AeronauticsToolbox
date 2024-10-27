% tri_score.fcn calculates the triaviation score given max speed, min speed
% and rate of climb
%
% INPUTS:
%   V_max: maximum speed (mph)
%   V_min: minimum (stall) speed (mph)
%   V_c: rate of climb (fpm)
%
% OUTPUTS:
%   score: triaviation score
%
% Written by:
%   Sam Jaeger
%   1/31/2023
%   Updated: 2/8/2023

function score = tri_score(V_max,V_min,V_c)
    score = (28110625*(V_max^2)*(V_c^2))/((4100625+(V_min^4))*10^9);
end