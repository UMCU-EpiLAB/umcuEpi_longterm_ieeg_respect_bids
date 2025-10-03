function [x,y,z,point_value] = improfile3D(image,p1,p2)

%--------------------------------------------------------------------------
%   
%   Returns coordinates and values of each voxel along a line between 
%   two points in a 3D image. 
%
%   INPUTS --
%   image : original image
%   p1 & p2 : endpoints of the line
%   p1 = [x1 y1 z1]
%   p2 = [x2 y2 z2]
%   OUTPUTS --
%   x, y, z : Coordinates of each voxels along the line
%   point_value : Value of each voxels [x y z]
%
% Author: Aurélien Sibellas
% e-mail: aurelien.sibellas@gmail.com
% Release: 1.0
% Release date: 18/03/2021
%--------------------------------------------------------------------------

% Euclidian distance
dist_euc = norm(p1 - p2);

% Number of intervalles
n_intervalle = round(dist_euc*2);
step = (p1 - p2)/n_intervalle;
pix_coords = zeros(n_intervalle, 3);

for cp = 0:n_intervalle
    pix_coords(cp+1, :) = round(p2 + cp*step);
end

point_list = sub2ind(size(image),pix_coords(:,1),pix_coords(:,2),pix_coords(:,3));
point_list = unique(point_list);

if isnan(point_list)
    point_value = [];
else
    point_value = image(point_list);
end

% Point coordinates
x = pix_coords(:,1);
y = pix_coords(:,2);
z = pix_coords(:,3);
