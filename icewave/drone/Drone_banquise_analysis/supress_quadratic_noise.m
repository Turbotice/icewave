function [V_s] = supress_quadratic_noise(V,x,y)

% This function supress the quadratic noise from a velocity field
% It returns the filtered velocity field
% It takes as arguments :
% - V : the velocity field [nx,ny,nt]
% - x : array of pixels along 1st dimension 
% - y : array of pixels along 2nd dimension

[nx,ny,n] = size(V);

[Y_quad,X_quad] = meshgrid(y,x);
V_s = zeros(nx,ny,n);
for i=1:n
    
    Vx = V(:,:,i);
    % fit by a quadratic function
    P = polyFit2D(Vx,X_quad,Y_quad,2,2);
    Pth = polyVal2D(P,X_quad,Y_quad,2,2);

    V_s(:,:,i) = Vx-Pth; % get rid off quaratic noise (drone movements)

end 