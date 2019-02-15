function [K,H] = surfaceCurvature(S,N)
%SURFACECURVATURE calculates the local curvature of each face in the input triangular mesh 
% 
% K = surfaceCurvature(surface,normals)
% 
% [K,H] = surfaceCurvature(surface,normals)
%
% This function will calculate an approximate curvature value for each face
% in the input triangular mesh. The surface should follow the format used
% by patch, where the surface is contained in a structure with the fields
% "faces" and "veritces"
% 
% Input: 
% 
%   surface - The surface to calculate curvature on, using the FV format
%   (Faces/vertices) used by patch, isosurface etc.
% 
%   normals - The normals of the surface at each vertex. These normals need
%   not be of unit length.
% 
%   Example:
%   To calculate the local curvature of an isosurface of an image, use the
%   following commands:
% 
%       s = isosurface(image,isoValue);
%
%       n = isonormals(image,s.vertices); 
% 
%       c = surfaceCurvature(s,n);
% 
%   Which can then be visualized with the command:
% 
%       patch(s,'FaceColor','flat','EdgeColor','none','FaceVertexCData',c)    
% 
% 
% Output:
% 
%   K = An Mx1 vector, where M is the number of faces, of the approximate
%   gaussian curvature at each face.
% 
%   H = An Mx1 vector, where M is the number of faces, of the approximate
%   mean curvature at each face.
% 
%
% References:
%
% [1] Theisel et al, "Normal Based Estimation of the Curvature Tensor for
% Triangular Meshes", Proceeeding of the Computer Graphics and
% Applications, 12th pacific Conference (PG '04)
% 
%Hunter Elliott 
%3/2010
%

if nargin < 2 || isempty(S) || isempty(N)
    error('Must input surface mesh and surface normals!')
end

if ~isfield(S,'vertices') || ~isfield(S,'faces')
    error('The input surface must be a structure using the FV format, with a field named vertices and a field named faces!')
end

%Number of faces
nTri = size(S.faces,1);

%Barycentric coordinates for location of interpolated normal
abc = ones(1,3) * 1/3; %This will estimate curvature at the barycenter of each face.

%Init array for curvature values
K = zeros(nTri,1);
H = zeros(nTri,1);

%Should probably vectorize/arrayfun this at some point...
for i = 1:nTri
    
    %Get the coordinates of this triangle's vertices
    X = S.vertices(S.faces(i,:),:);
    
    %Get the normal vectors for these vertices
    n = N(S.faces(i,:),:);
    
    %Interpolated normal
    ni = abc * n;
    
    %Triangle normal
    m = cross(X(2,:)-X(1,:),X(3,:)-X(2,:));
    
    %Gaussian curvature - formula (12) in ref [1]
    K(i) = det(n) / (dot(ni,ni)*dot(ni,m));
    
           
    %H from formula (13) in ref [1]
    h = cross(n(1,:),X(3,:)-X(2,:))+...
        cross(n(2,:),X(1,:)-X(3,:))+...
        cross(n(3,:),X(2,:)-X(1,:));
    
    %Mean curvature - formula 
    H(i) = .5*dot(ni,h) / (sqrt(dot(ni,ni))*dot(ni,m));
              
end

