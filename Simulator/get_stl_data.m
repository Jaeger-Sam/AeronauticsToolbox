% get_stl_data1.fcn opens an stl file and returns a matrix of verticies and
% normal vectors of each panel. Adapated from Maziar Hemati NAflightsim
% code. Assumes that input stl file is in meters.
%
% [N,v,n] = get_stl_data(filen,units)
%
% Input:
%   filen: [string] file name
%   units: if units == 'US' then will convert verticies to feet. Otherwise
%       
% Output:
%   N: number of panels
%   v: [Nx3] matrix of vertices (v1,v2,v3)_i for i=1,...,N
%   n: [Nx3] matrix of normal vectors (n1,n2,n3)_i for i=1,...,N
%
% Sam Jaeger, Maziar Hemati
% 1/30/2024
%   Revised: 2/27/2024
%   Revised: 10/15/2024

function [N,v,n] = get_stl_data(filen,units)
    file_id = fopen(filen);
    C=textscan(file_id,'%s');
    kk=cellfun(@length,C);
    
    % count number of faces
    %   each face should only have 1 normal
    N = 0;
    for ii = 1:kk
        if strcmp('normal',C{1}{ii})
            N = N+1;
        end
    end
    
    n=NaN(N,3); v=NaN(N,9);
    
    count = 0;
    vcount = 0;
    for ii = 1:kk
        if strcmp('normal',C{1}{ii})
            count = count+1;
            for jj = 1:3
                n(count,jj) = str2double(C{1}{ii+jj});
            end
        end
    
        if strcmp('vertex',C{1}{ii})
            for jj = 1:3
                v(count,3*vcount+jj) = str2double(C{1}{ii+jj});
            end
            vcount= vcount+1;
            if vcount==3
                vcount=0;
            end
        end
    
    end
    
    if units == 'US'
        v = v*3.28084;
    end
    
    save([filen '.mat'],'N','n','v')
end