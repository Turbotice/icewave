
clear all
close all

%% Load data
base = 'W:/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/matData/DJI_0402/Data_ice_Dt4_W64/';

filename = 'PIV_processed_Dt4_b1_DJI_0402_post_processed.mat' ;

disp('Loading data...')
load([base filename])
disp('Data loaded')

%% Define fig_folder
fig_folder = 'W:/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/Figures_report/DJI_0402/';

if exist(fig_folder,'dir') ~= 7
    mkdir(fig_folder)
end

%% Scales

facq_x = 1/m.fx; % scale in box / meter
facq_t = 1/m.ft; % scale in frame / sec
scale_V = m.scale_V; % factor scaling for velocity in meter / s
W = m.w ;
font_size = 13;

%% Plot raw velocity field 
raw_figure = figure(1);

frame_idx = 5000;
V_raw = flip(m.Vx(1:105,:,frame_idx)*scale_V,1);
x = (1:1:size(V_raw,1));
y = (1:1:size(V_raw,2));
pcolor(x/facq_x,y/facq_x,V_raw')
shading interp
%     colormap(redblue)
xlabel('$x$ (m)');
ylabel('$y$ (m)');
shading interp
% axis([0 size(V_raw,2)/fx 0 size(X,1)/fx])
axis image
%     set_Papermode(gcf)
cbar = colorbar();
cbar.Label.String = '$ V_x(x,y) \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'Latex';
%     cbar.Label.FontSize = font_size;

% caxis([-caxis_amp caxis_amp])
ax = gca;
ax.FontSize = 13;
set_Papermode(gcf)

filename = ['Vx_raw_velocity_field_frame_' num2str(frame_idx)];
saveas(raw_figure,[fig_folder filename],'fig')
saveas(raw_figure,[fig_folder filename],'pdf')
saveas(raw_figure,[fig_folder filename],'png')

%% Supress Drone motion 
[nx,ny,nt] = size(m.Vx);
x = (1:1:nx);
y = (ny:-1:1);

% compute the mean component of the velocity field for each frame
Vxmoy = mean(mean(m.Vx,2),1);
Vymoy = mean(mean(m.Vy,2),1);
% reduce the velocity field from its mean value
% m.Vx = m.Vx - Vxmoy;
% m.Vy = m.Vy - Vymoy;

Vx = supress_quadratic_noise(m.Vx*scale_V,x,y);
Vy = supress_quadratic_noise(m.Vy*scale_V,x,y);

disp('Drone motion corrected')

%% Compute time Fourier transform of raw Vx

disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;
[FFT_t_raw,TF_spectrum_raw,f_raw] = temporal_FFT(m.Vx(:,:,:)*scale_V,padding_bool,add_pow2,facq_t);

%%
fig_spectrum = figure(1); 
loglog(f_raw,TF_spectrum_raw)
grid on 
xlabel('$f \: \rm (Hz)$','Interpreter','latex')
ylabel('$\overline {\langle V_x \rangle _{x,y}} (f) \: \rm (m.s^{-1})$','Interpreter','latex')
axis([1e-3 30 1e-4 6e-2])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;

file = [fig_folder 'FFT_spectrum_raw_velocity_field'];
saveas(fig_spectrum,file,'fig')
saveas(fig_spectrum,file,'pdf')
saveas(fig_spectrum,file,'png')
%% Compute time Fourier transform of raw Vx

disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;
[FFT_t,TF_spectrum,f] = temporal_FFT(Vx,padding_bool,add_pow2,facq_t);

%%
fig_spectrum = figure(2); 
loglog(f_raw,TF_spectrum_raw)
hold on 
loglog(f,TF_spectrum)
grid on
xlabel('$f \: \rm (Hz)$','Interpreter','latex')
ylabel('$\overline {\langle V_x \rangle _{x,y}} (f) \: \rm (m.s^{-1})$','Interpreter','latex')
legend('no corrections','quadratic corrections','Location','southwest')
axis([1e-3 30 1e-4 6e-2])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;

file = [fig_folder 'FFT_spectrum_velocity_field_superposition_quadratic_correction'];
saveas(fig_spectrum,file,'fig')
saveas(fig_spectrum,file,'pdf')
saveas(fig_spectrum,file,'png')

%% Plot Demodulated fields 

disp('Getting demodulated fields')
selected_freq = [0.2 0.7];
x_bound = [1 105];
caxis_amp = 8e-3; % amplitude of the colorbar in meter/second
fig_name = ['Demodulated_field_Vx_caxis_' num2str(caxis_amp) 'ms'];

save_image = 1;
save_video = 0;
plot_demodulated_field(FFT_t,f,facq_x,selected_freq,x_bound,caxis_amp,1,fig_folder,fig_name,save_image,save_video)

%% Get attenuation coefficient 

% indices used to fit the exponential decay
i_min = 1;
i_max = 105;

disp('Getting attenuation coefficient')
selected_freq = [0.1 0.7]; % selected frequencies between which we proceed to the analysis
x_bound = [12 67]; % selected boundaries at which we perform 2D FFT
freq_thresh = [0.46];
x_bound_thresh = round([27 67]*facq_x);

new_folder_fig = [fig_folder 'Attenuation_fit/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

fig_name = 'Decay_law_video_mode_1';
save_image = 1;
save_video = 0;
left_bool = 0;

[lambda,dist_fit,freq] = get_attenuation_coef(FFT_t,f,facq_x,selected_freq,x_bound,freq_thresh,x_bound_thresh,left_bool,new_folder_fig,fig_name,save_image,save_video);

% save data to plot attenuation coefficient
freq_min_txt = replace(num2str(selected_freq(1)),'.','p');
freq_max_txt = replace(num2str(selected_freq(2)),'.','p');
f_thresh_txt = replace(num2str(freq_thresh(1)),'.','p');
attenuation_file = [fig_folder 'attenuation_coef_data_fmin' freq_min_txt '_fmax' freq_max_txt '_fthresh' f_thresh_txt];
%attenuation_file = replace(attenuation_file,'.','p');
attenuation_file = [attenuation_file '_mode_1'];
save(attenuation_file,'lambda','dist_fit','freq','selected_freq','x_bound','freq_thresh','x_bound_thresh')

disp('DONE.')