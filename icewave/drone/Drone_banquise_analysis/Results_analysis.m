%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

%base = 'W:/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/matData/DJI_0402/Data_ice_Dt4_W64/';
base = 'E:/Rimouski_2024/Data/2024/0219/matData/waves_012/';

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


%%
% %% Scaling 
% fe = 29.97; % Frame rate in Hz
% nb_pass = m.s_param{6,2};
% W = 32; % size of the window for DIC algorithm
% font_size = 13;
% 
% % ##########################################
% L_x = 3840; % size of the image in pixel, along x-axis
% h_drone = 330*0.3048; % height of the drone in meter
% theta_x = 32.75; % semi AFOV of the drone, along x-axis, in °
% 
% fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
% fx = fx_pix*2/W;
% % fx = 1;
% % ##########################################
% 
% % fx = 0.8857; % for W = 64; 
% %fx = 0.8857*2; % spatial scale in boxes/meter
% Dt = 2; % step between two frames that were compared during the PIV algorithm 
% 
% scale_V = (fe/Dt) / fx_pix; % scale of the velocity in m/s
% scale_V = 1;

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


Vx = supress_quadratic_noise(m.Vx,x,y);
Vy = supress_quadratic_noise(m.Vy,x,y);

disp('Drone motion corrected')
% m_corrected = m;
% m_corrected.Vx = Vx;
% m_corrected.Vy = Vy;

% save structure with corrected velocity
% scale_file_name = replace(matname,'_total_processed.mat','_corrected.mat');
% save(scale_file_name,'m_corrected','-v7.3');

%% Get velocity map 
[nx,ny,nt] = size(m.Vx);

% create a meshgrid
y = (ny:-1:1);
x = (1:1:nx);

[X,Y]=meshgrid(x,y);

i = 5000;
caxis_amp = 1.0;
x_max = nx; % max value for ice 
facq_x = 1/m.fx;

figure,
surf(X./facq_x,Y./facq_x,m.Vx(1:x_max,:,i)')
view(2)
shading interp
axis image
axis([0 size(X,2)/facq_x 0 size(Y,1)/facq_x]);
caxis([-caxis_amp caxis_amp])
cbar =colorbar();
cbar.Label.String = '$V_x \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
cbar.FontSize = 13;
%title('$V_x$','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
ax = gca;
ax.FontSize = 13;

%% Plot main features of the velocity field
disp('Get scaled velocity')
idx_frame = 20;
caxis_amp = 1.0; % Amplitude of the colorbar scale (in meter / second)

plot_velocity_features(Vx,Vy,facq_x,idx_frame,caxis_amp,fig_folder);

%% Get a movie of the velocity field 
caxis_amp = [-1 1]; %% modifier avec le quantile de 0.9
fig_name = [fig_folder 'Velocity_field_Vx_movie'];
fps = facq_t ;
movie_velocity_field(Vx,facq_x,facq_t,caxis_amp,fps,fig_name)

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
selected_freq = [0.07 1.0];
x_bound = [1 size(FFT_t,1)];
caxis_amp = 0.04; % amplitude of the colorbar in meter/second
fig_name = 'Demodulated_field_continuous';

save_image = 0;
save_video = 1;
plot_demodulated_field(FFT_t,f,facq_x,selected_freq,x_bound,caxis_amp,fig_folder,fig_name,save_image,save_video)

%% Get wave vectors 
disp('Getting wave vectors')
selected_freq = [0.1 1.0]; % selected frequencies between which we proceed to the analysis
x_bound = [1 48*facq_x]; % selected boundaries at which we perform 2D FFT
padding_bool = 1;
add_pow2 = 2; % additional power of 2 for padding 
black_mask = 10; % in pixels
caxis_amp = -1;
fig_name = 'Spatial_Fourier_space_video';

save_image = 0;
save_video = 1;


[k,freq] = get_wave_vectors(FFT_t,f,facq_x,selected_freq,x_bound,padding_bool,add_pow2,black_mask,caxis_amp,fig_folder,fig_name,save_image,save_video);

% save data to plot dispersion relation
filename = ['dispersion_relation_data_continuous_fmin' num2str(selected_freq(1)) '_fmax' num2str(selected_freq(2)) '_add_pow' num2str(add_pow2)];
filename = replace(filename,'.','p');
dispersion_file = [fig_folder filename];
save(dispersion_file,'k','freq','selected_freq','x_bound','add_pow2','black_mask')

%% Get attenuation coefficient 
disp('Getting attenuation coefficient')
selected_freq = [0.07 1.1]; % selected frequencies between which we proceed to the analysis
x_bound = [round(3*facq_x) size(Vx,1)]; % selected boundaries at which we perform 2D FFT
freq_thresh = 0.8;
x_bound_thresh = [round(3*facq_x) round(40*facq_x)];

new_folder_fig = [fig_folder 'Attenuation_fit/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

fig_name = 'Decay_law_video';
save_image = 0;
save_video = 1;

[lambda,dist_fit,freq] = get_attenuation_coef(FFT_t,f,facq_x,selected_freq,x_bound,freq_thresh,x_bound_thresh,new_folder_fig,fig_name,save_image,save_video);

% save data to plot attenuation coefficient
freq_min_txt = replace(num2str(selected_freq(1)),'.','p');
freq_max_txt = replace(num2str(selected_freq(2)),'.','p');
f_thresh_txt = replace(num2str(freq_thresh),'.','p');
attenuation_file = [base_fig 'attenuation_coef_data_fmin' freq_min_txt '_fmax' freq_max_txt '_fthresh' f_thresh_txt];
%attenuation_file = replace(attenuation_file,'.','p');
save(attenuation_file,'lambda','dist_fit','freq','selected_freq','x_bound','freq_thresh','x_bound_thresh')

disp('DONE.')