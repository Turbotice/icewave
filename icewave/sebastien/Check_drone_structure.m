
clear all

%% Check structure of drone data 

base = 'W:/SagWin2024/Data/0211/Drones/bernache/matData/18-stereo_001/Data_not_structured/';
filename = 'Data_PIV_oblique_Dt4_W32_bernache.mat';

disp('Loading data..')
load([base filename])
disp('Data loaded')

%% Test backward_projection function to get vertical velocity 
[nx,ny,nt] = size(s.Vx);
x = (1:nx);
y = (1:ny);
t = (1:nt);

W = s.param.W;
x_pix = ((W + 1)/2 : W/2 : (nx*W + 1)/2); % x-index in pixel system of each box 
y_pix = ((W + 1)/2 : W/2 : (ny*W + 1)/2); % y-index in pixel system of each box

Fy = griddedInterpolant({x_pix,y_pix,t},s.Vy); % in pixels

%%
y_0 = (2160 + 1)/2;

[Fz] = backward_projection(Fy,y_0,s.param.H,s.param.alpha_0,s.param.focale,s.param.fps,s.param.Dt);

%%
[X_pix,Y_pix,T] = meshgrid(x_pix,y_pix,t);
Vz = Fz(X_pix,Y_pix,T);

%% 
Vz = permute(Vz,[2,1,3]);
Vz = - Vz;

[X_pix,Y_pix] = meshgrid(x_pix,y_pix);
[Xreal,Yreal] = projection_real_space();