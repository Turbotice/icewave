%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

base = '//192.168.1.70/Share/Data/0215/Drones/mesange/matData/waves_003/';
%base = 'E:/PIVlab_drone/matdata/DJI_0308_Dt4_W64_full/';

filename = 'PIV_processed_i01_Dt2_b1_W32_full_total_processed.mat';
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

% %% Scaling 
% fe = 29.97; % Frame rate in Hz
% nb_pass = m.s_param{6,2};
% W = 36; % size of the window for DIC algorithm
% font_size = 13;
% 
% % ##########################################
% L_x = 3840; % size of the image in pixel, along x-axis
% h_drone = 140; % height of the drone in meter
% theta_x = 32.75; % semi AFOV of the drone, along x-axis, in °
% 
% fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
% fx = fx_pix*2/W;
% % ##########################################
% 
% % fx = 0.8857; % for W = 64; 
% %fx = 0.8857*2; % spatial scale in boxes/meter
% Dt = 4; % step between two frames that were compared during the PIV algorithm 
% 
% scale_V = (fe/Dt) / fx_pix; % scale of the velocity in m/s
% % store scales in structure m
% m.scale_V = scale_V;
% m.ft = fe;
% m.fx = fx;
%% Get Histogram displacement 

% Vx = m.Vx(10:end,:,:) - mean(mean(m.Vx(10:end,:,:),2),1);
filename = [fig_folder 'histogramm_displacements'];

get_histogram_displacement(m.Vx/scale_V,W,font_size,filename);

%% Get profile of velocity
i_x = 20;
i_y = 20;
plot_located_profile(m.Vx,i_x,i_y,facq_x,facq_t,fig_folder);

disp(mean(abs(m.Vx(i_x,i_y,:))))


%% Get rid off quadratic noise (drone movements)
[nx,ny,nt] = size(m.Vx);
x = (1:1:nx);
y = (ny:-1:1);

% compute the mean component of the velocity field for each frame
Vxmoy = mean(mean(m.Vx,2),1);
Vymoy = mean(mean(m.Vy,2),1);
% reduce the velocity field from its mean value
% m.Vx = m.Vx - Vxmoy;
% m.Vy = m.Vy - Vymoy;

m_corrected = m;
Vx = supress_quadratic_noise(m.Vx,x,y);
Vy = supress_quadratic_noise(m.Vy,x,y);

m_corrected.Vx = Vx;
m_corrected.Vy = Vy;

% save structure with corrected velocity
scale_file_name = replace(matname,'.mat','_corrected.mat');
save(scale_file_name,'m_corrected','-v7.3');

%% Plot main features of the velocity field
disp('Get scaled velocity')
idx_frame = 20;
caxis_amp = 1.0; % Amplitude of the colorbar scale (in meter / second)

plot_velocity_features(Vx,Vy,facq_x,scale_V,idx_frame,caxis_amp,fig_folder);

%% Get time Fourier transform
disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;
[FFT_t,TF_spectrum,f] = temporal_FFT(Vx(:,:,:),padding_bool,add_pow2,facq_t);

%%
fig_spectrum = figure; 
loglog(f,TF_spectrum)
grid on 
xlabel('$f \: \rm (Hz)$','Interpreter','latex')
ylabel('$\overline {\langle V_x \rangle _{x,y}} (f) \: \rm (m.s^{-1})$','Interpreter','latex')
ax = gca;
ax.FontSize = 13;

file = [fig_folder 'FFT_spectrum'];
saveas(fig_spectrum,file,'fig')

%% Get demodulated field
disp('Getting demodulated fields')
selected_freq = [0.05 0.7];
x_bound = [1 50*facq_x];
caxis_amp = -1; % amplitude of the colorbar in meter/second
fig_name = 'Demodulated_field_continuous';

save_image = 0;
save_video = 1;
plot_demodulated_field(FFT_t,f,facq_x,selected_freq,x_bound,caxis_amp,fig_folder,fig_name,save_image,save_video)

%% Get wave vectors 
disp('Getting wave vectors')
selected_freq = [0.05 0.7]; % selected frequencies between which we proceed to the analysis
x_bound = [1 50*facq_x]; % selected boundaries at which we perform 2D FFT
padding_bool = 1;
add_pow2 = 2; % additional power of 2 for padding 
black_mask = 10;
caxis_amp = -1;
fig_name = 'Spatial_Fourier_space_continuous_video';

save_image = 0;
save_video = 1;


[k,freq] = get_wave_vectors(FFT_t,f,facq_x,selected_freq,x_bound,padding_bool,add_pow2,black_mask,caxis_amp,fig_folder,fig_name,save_image,save_video);

% save data to plot dispersion relation
filename = ['dispersion_relation_data_continuous_fmin' num2str(selected_freq(1)) '_fmax' num2str(selected_freq(2)) '_add_pow' num2str(add_pow2)];
filename = replace(filename,'.','p');
dispersion_file = [base_fig filename];
save(dispersion_file,'k','freq','selected_freq','x_bound','add_pow2','black_mask')

%% Get attenuation coefficient 
disp('Getting attenuation coefficient')
selected_freq = [0.1 0.7]; % selected frequencies between which we proceed to the analysis
x_bound = [1 size(Vx(:,:,:),1)]; % selected boundaries at which we perform 2D FFT
freq_thresh = 0.35;
x_bound_thresh = [1 40];

new_folder_fig = [fig_folder 'Attenuation_fit/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

fig_name = 'Decay_law_video';
save_image = 1;
save_video = 1;

[lambda,dist_fit,freq] = get_attenuation_coef(FFT_t,f,facq_x,selected_freq,x_bound,freq_thresh,x_bound_thresh,new_folder_fig,fig_name,save_image,save_video);

% save data to plot attenuation coefficient

attenuation_file = [base_fig 'attenuation_coef_data_fmin' num2str(selected_freq(1)) '_fmax' num2str(selected_freq(2)) '_f_thresh' num2str(freq_thresh)];
attenuation_file = replace(attenuation_file,'.','p');
save(attenuation_file,'lambda','dist_fit','freq','selected_freq','x_bound','freq_thresh','x_bound_thresh')

disp('DONE.')