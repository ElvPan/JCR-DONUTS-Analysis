function [e_mask]=ellipsoid2(ix,iy,iz,cx,cy,cz,a,b,c)

% Draws a ellipsoid mask with centre cx,cy,cz in image ix,iy,iz with axes
% a, b,c

[x,y,z]=meshgrid(-(ix/2):(ix/2-1),-(iy/2):(iy/2-1),-(iz/2):(iz/2-1));
e_mask=((((x-cx)/a).^2+((y-cy)/b).^2+((z-cz)/c).^2)<=1);