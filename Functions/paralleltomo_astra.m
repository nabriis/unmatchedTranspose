function [A,A_m] = paralleltomo_astra(N,theta,p,d,proj_type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Object geometry
vol_geom  = astra_create_vol_geom(N,N); %astra_create_vol_geom(N,N,-domain/2,domain/2,-domain/2,domain/2);

% Projection geometry
proj_geom = astra_create_proj_geom('parallel', d/p, p, degtorad(theta));

% Create the Spot operator for ASTRA using the GPU.
P = opPermutation(reshape(reshape(1:p*length(theta),p,length(theta))',p*length(theta),1));
A = opFoG(P',opTomo(proj_type, proj_geom, vol_geom));

if nargout > 1
    angles_astra = -(theta + 90)*pi/180;
    proj_geom = astra_create_proj_geom('parallel', d/p, p, angles_astra);
    proj_id1 = astra_create_projector(proj_type, proj_geom, vol_geom);
    
    % Generate the projection matrix for this projection model.
    % This creates a matrix W where entry w_{i,j} corresponds to the
    % contribution of volume element j to detector element i.
    matrix_id1 = astra_mex_projector('matrix', proj_id1);
    %
    % % Get the projection matrix as a Matlab sparse matrix.
    A_m = astra_mex_matrix('get', matrix_id1);
end

end

