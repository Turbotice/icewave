%% This script is used to filter and analyze data from the 02/11 Saguenay, with Oblique view

clear all
close all

%% Load Data 

base = 'W:/SagWin2024/Data/0211/Drones/bernache/matData/18-stereo_001/';

disp('Loading data..')
load([base 'Data_PIV_oblique_Dt4_W32_bernache.mat'])
disp('Data loaded')

%% Try to apply a median filter to the data (median in space)
[nx,ny,nt] = size(s.dz);

box_filter = [5 5]; % number of neighbours used to apply the median filter

filter_dz = zeros(nx,ny,nt);
for i0 = 1 : nt
    disp(i0)
    % Apply a median filter 
    img = s.dz(:,:,i0);
    J = medfilt2(img,box_filter);
    filter_dz(:,:,i0) = J;
    
end 


%% Try to apply a median filter in  both space and time 
[nx,ny,nt] = size(s.dz);

box_filter = [3 3 3]; % number of neighbours used to apply the median filter

filter_dz = medfilt3(s.dz,box_filter);
disp('Median filter applied')
%%

%# Buoy #1
% obj_txt = 'buoy_1';
% x_bernache = 1284;
% y_bernache = 752;

%# Buoy #5
obj_txt = 'buoy_5';
x_bernache = 3204;
y_bernache = 759;

x_mesange = 586;
y_mesange = 911;

signal_bernache = get_1D_signal(x_bernache,y_bernache,s.dz,32);
filter_sig = get_1D_signal(x_bernache,y_bernache,filter_dz,32);

figure, 
plot(s.t,signal_bernache)
hold on 
plot(s.t,filter_sig)
grid on 

%% 
i = 950;
figure, 
surf(s.X,s.Y,s.dz(:,:,i))
shading interp
view(2)

figure, 
surf(s.X,s.Y,filter_dz(:,:,i))
shading interp
view(2)


%% Superposition of drone data and buoys 

%% ############# FUNCTION SECTION ###############

function [signal] = get_1D_signal(x_obj,y_obj,dz,W)

% This function computes the vertical velocity at a given position of the
% PIV field, and returns it as a 1D signal
% The function takes the following arguments : 
% - x_obj : x-coordinate on the camera (in pixels)
% - y_obj : y-coordinate on the camera (in pixels)
% - dz : vertical velocity matrix [nx,ny,nt]
% - W : window size used for the PIV processing (pixels)

% x_0 = 3840/2;
% y_0 = 2160/2;

% object position in the sea ice framework
% Y_b = (y_obj - y_0)*H/sin(alpha_0)./(f*sin(alpha_0) + (y_obj - y_0).*cos(alpha_0));
% X_b = (x_obj - x_0)*H./(f*sin(alpha_0) + (y_obj - y_0).*cos(alpha_0));
% 
% Y_b = - Y_b;

% get closest point in the real framework 
idx_x = floor(x_obj*2/W);
idx_y = floor(y_obj*2/W);

signal = squeeze(dz(idx_x,idx_y,:));


end