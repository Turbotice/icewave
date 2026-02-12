%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

base = 'W:/SagWin2024/Data/0226/Drones/mesange/matData/frac_001/';

% base = 'E:/Rimouski_2024/Data/2024/0219/matData/waves_012/';

filename = 'PIV_processed_i00_Dt3_b1_W32_full_total_processed.mat';
matname = [base filename];
%%
disp('Loading Data..');
load(matname);
disp('Data loaded');
%% Creates a folder where to save the generated plots
base_fig = base;
%base_fig = 'E:/PIVlab_drone/';
fig_folder = [base_fig 'Plots/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Scales

facq_x = 1/m.fx; % scale in box / meter
facq_t = 1/m.ft; % scale in frame / sec
scale_V = m.scale_V; % factor scaling for velocity in meter / s
W = m.w ;
font_size = 13;

%% Scaling 
% fe = 29.97; % Frame rate in Hz
% nb_pass = m.s_param{6,2};
% W = 32; % size of the window for DIC algorithm
% font_size = 13;
% 
% % ##########################################
% L_x = 3840; % size of the image in pixel, along x-axis
% h_drone = 100.3; % height of the drone in meter
% theta_x = 34.15; % semi AFOV of the drone, along x-axis, in °
% 
% fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
% fx = fx_pix*2/W;
% % fx = 1;
% % ##########################################
% 
% % fx = 0.8857; % for W = 64; 
% %fx = 0.8857*2; % spatial scale in boxes/meter
% Dt = 4; % step between two frames that were compared during the PIV algorithm 
% 
% scale_V = (fe/Dt) / fx_pix; % scale of the velocity in m/s


%% store scales in structure m
% m.scale_V = scale_V;
% m.ft = 1/fe;
% m.fx = 1/fx;
% m.h_drone = h_drone;
% m.theta = theta_x;
%%
m = m_scaled;
%% Get Histogram displacement 

% Vx = m.Vx(10:end,:,:) - mean(mean(m.Vx(10:end,:,:),2),1);
filename = [fig_folder 'histogramm_displacements_Vx_time_average'];
average_bool = 1;
get_histogram_displacement(m.Vx/scale_V,W,average_bool,font_size,filename);

%% Get profile of velocity
i_x = 150;
i_y = 20;
plot_located_profile(m.Vx,i_x,i_y,facq_x,facq_t,fig_folder);

disp(mean(abs(m.Vx(i_x,i_y,:))))
