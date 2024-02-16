%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

base = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/Traitement_drone_20230310/matData/test/';
%base = 'E:/PIVlab_drone/matdata/DJI_0308_Dt4_W64_full/';

filename = 'PIV_processed_Dt4_b1_W64_full_test_images_total_processed.mat';
matname = [base filename];
%%
disp('Loading Data..');
load(matname);
disp('Data loaded');
%% Creates a folder where to save the generated plots
base_fig = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/';
%base_fig = 'E:/PIVlab_drone/';
fig_folder = [base_fig 'DJI_0308_figures_Dt4_W64_full/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Scaling 
fe = 29.97; % Frame rate in Hz
nb_pass = m.s_param{6,2};
W = 2^(9-nb_pass); % size of the window for DIC algorithm
font_size = 13;

% ##########################################
L_x = 3840; % size of the image in pixel, along x-axis
h_drone = 78.03; % height of the drone in meter
theta_x = 32.75; % semi AFOV of the drone, along x-axis, in °

fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
fx = fx_pix*2/W;
% ##########################################

% fx = 0.8857; % for W = 64; 
%fx = 0.8857*2; % spatial scale in boxes/meter
Dt = 4; % step between two frames that were compared during the PIV algorithm 

scale_V = (fe/Dt) / fx_pix; % scale of the velocity in m/s
% store scales in structure m
m.scale_V = scale_V;
m.ft = fe;
m.fx = fx;
%% Get Histogram displacement 

% Vx = m.Vx(10:end,:,:) - mean(mean(m.Vx(10:end,:,:),2),1);
get_histogram_displacement(m.Vx,W,fig_folder,font_size);

%% Get profile of velocity
i_x = 20;
i_y = 20;
plot_located_profile(m.Vx,i_x,i_y,fx,fe,fig_folder);

disp(mean(abs(m.Vx(i_x,i_y,:))))


%% Get rid off quadratic noise (drone movements)
[nx,ny,nt] = size(m.Vx);
x = (1:1:nx);
y = (ny:-1:1);

% compute the mean component of the velocity field for each frame
Vxmoy = mean(mean(m.Vx,2),1);
Vymoy = mean(mean(m.Vy,2),1);
% reduce the velocity field from its mean value
m.Vx = m.Vx - Vxmoy;
m.Vy = m.Vy - Vymoy;

m.Vx = supress_quadratic_noise(m.Vx,x,y);
m.Vy = supress_quadratic_noise(m.Vy,x,y);

%% Get Velocity scaled
disp('Get scaled velocity')
idx_frame = 1;
caxis_amp = 1.0; % Amplitude of the colorbar scale (in meter / second)
[Vx_s,Vy_s] = get_scaled_velocity_field(m.Vx,m.Vy,fx,scale_V,idx_frame,caxis_amp,fig_folder);

m_scaled = m;
m_scaled.Vx = Vx_s;
m_scaled.vy = Vy_s;

% save structure with scaled velocity
scale_file_name = replace(matname,'.mat','_scaled.mat');
save(scale_file_name,"m_scaled");

%% Get time Fourier transform
disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;
[FFT_t,TF_spectrum,f] = temporal_FFT(m_scaled.Vx(10:end,:,:),padding_bool,add_pow2,fe);

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
selected_freq = [0.1 0.8];
x_bound = [1 size(m.Vx(10:end,:,:),1)];
caxis_amp = -1; % amplitude of the colorbar in meter/second
plot_demodulated_field(FFT_t,f,fx,selected_freq,x_bound,caxis_amp,fig_folder)

%% Get wave vectors 
disp('Getting wave vectors')
selected_freq = [0.1 0.8]; % selected frequencies between which we proceed to the analysis
x_bound = [1 size(m.Vx(10:end,:,:),1)]; % selected boundaries at which we perform 2D FFT
padding_bool = 1;
add_pow2 = 2; % additional power of 2 for padding 
black_mask = 10;
caxis_amp = -1;

[k,freq] = get_wave_vectors(FFT_t,f,m.fx,selected_freq,x_bound,padding_bool,add_pow2,black_mask,caxis_amp,fig_folder);

% save data to plot dispersion relation
dispersion_file = [base_fig 'dispersion_relation_data_fmin' num2str(selected_freq(1)) '_fmax' num2str(selected_freq(2)) '_add_pow' num2str(add_pow2)];
dispersion_file = replace(dispersion_file,'.','p');
save(dispersion_file,'k','freq','selected_freq','x_bound','add_pow2','black_mask')

%% Get attenuation coefficient 
disp('Getting attenuation coefficient')
selected_freq = [0.1 0.7]; % selected frequencies between which we proceed to the analysis
x_bound = [1 size(m.Vx(10:end,:,:),1)]; % selected boundaries at which we perform 2D FFT
freq_thresh = 0.35;
x_bound_thresh = [1 40];

new_folder_fig = [fig_folder 'Attenuation_fit/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

[lambda,dist_fit,freq] = get_attenuation_coef(FFT_t,f,fx,selected_freq,x_bound,freq_thresh,x_bound_thresh,new_folder_fig);

% save data to plot attenuation coefficient

attenuation_file = [base_fig 'attenuation_coef_data_fmin' num2str(selected_freq(1)) '_fmax' num2str(selected_freq(2)) '_f_thresh' num2str(freq_thresh)];
attenuation_file = replace(attenuation_file,'.','p');
save(attenuation_file,'lambda','dist_fit','freq','selected_freq','x_bound','freq_thresh','x_bound_thresh')

disp('DONE.')