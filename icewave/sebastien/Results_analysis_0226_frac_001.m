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
i_x = 30;
i_y = 60;
profile_fig = plot_located_profile(m.Vx,i_x,i_y,facq_x,facq_t,0);
set_Papermode(profile_fig)

disp(mean(abs(m.Vx(i_x,i_y,:))))

fig_file_name = [fig_folder 'Time_profile_ix_' num2str(i_x) '_iy_' num2str(i_y)];
saveas(profile_fig,fig_file_name,'fig')

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

%% Plot main features of the velocity field
disp('Get scaled velocity')
idx_frame = 20;
caxis_amp = 1.0; % Amplitude of the colorbar scale (in meter / second)

plot_velocity_features(Vx,Vy,facq_x,idx_frame,caxis_amp,fig_folder,0);

%% Get a movie of the velocity field 
% computes Quantile of the velocity field 
Q = quantile(m.Vx,[0.1 0.9],'all');

%%
fig_folder_E = 'E:/Rimouski_2024/Data/2024/0226/Drones/mesange/12-FRAC_001/';
caxis_amp = [Q(1)-0.2 Q(2)+0.2]; 
fig_name = [fig_folder_E 'Velocity_field_Vx_movie_test'];
fps = facq_t ;
step = 5; % prints every 5 frames 
movie_velocity_field(m.Vx,facq_x,facq_t,caxis_amp,fps,fig_name,0)

%% Show quickly the velocity field 
video_filename = [fig_name '.avi']; % folder where the video is saved
vid = VideoWriter(video_filename);
vid.FrameRate = 10;
open(vid)

figure, 
x = (nx:-1:1)*m.fx;
y = (ny:-1:1)*m.fx;
t = (1:1:nt)*m.ft;
relevant_idx = (1:3:nt);

for i0 = 1:length(relevant_idx)
    i = relevant_idx(i0);
    
    pcolor(x,y,Vx(:,:,i)');
    xlabel('$x \: \rm (m)$')
    ylabel('$y \: \rm (m)$')
    title(['$t = ' num2str(t(i)) ' \: \rm (s)$'])

    shading interp
    axis image
    colormap(redblue)
    colorbar()
    caxis([-0.6 0.6])
    set(gca,'XDir','reverse')
    T(i0)=getframe(gcf);
%     pause(0.1)
end 

writeVideo(vid,T)
close(vid)  

%% Compute Fourier Transform 
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
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;

file = [fig_folder 'FFT_spectrum'];
saveas(fig_spectrum,file,'fig')

%% Get demodulated field
disp('Getting demodulated fields')
selected_freq = [0.1 0.9];
x_bound = [1 size(FFT_t,1)];
caxis_amp = 0.01; % amplitude of the colorbar in meter/second
fig_name = 'Demodulated_field_Vx';

save_image = 0;
save_video = 1;
plot_demodulated_field(FFT_t,f,facq_x,selected_freq,x_bound,caxis_amp,0,fig_folder,fig_name,save_image,save_video)

%% Get wave vectors 
disp('Getting wave vectors')
selected_freq = [0.1 0.9]; % selected frequencies between which we proceed to the analysis
x_bound = [1 size(FFT_t,1)]; % selected boundaries at which we perform 2D FFT
padding_bool = 1;
add_pow2 = 2; % additional power of 2 for padding 
black_mask = 10; % in pixels
caxis_amp = -1;
fig_name = 'Spatial_Fourier_space_video';

save_image = 1;
save_video = 1;


[k,freq] = get_wave_vectors(FFT_t,f,facq_x,selected_freq,x_bound,padding_bool,add_pow2,black_mask,caxis_amp,fig_folder,fig_name,save_image,save_video);

% save data to plot dispersion relation
filename = ['dispersion_relation_data_Vx_fmin' num2str(selected_freq(1)) '_fmax' num2str(selected_freq(2)) '_add_pow' num2str(add_pow2)];
filename = replace(filename,'.','p');
dispersion_file = [fig_folder filename];
save(dispersion_file,'k','freq','selected_freq','x_bound','add_pow2','black_mask')

%% Quick plot of dispersion relation 

mask = (omega < 4.4) & (omega > 1.08);
filt_k = k(mask);
filt_omega = omega(mask);

g = 9.81;
k_list = linspace(0.07,5,100);
deep_water = sqrt(g*k_list);

omega = 2*pi*freq;
figure,
loglog(filt_k,filt_omega,'o')
grid on 
hold on 
loglog(k_list,deep_water,'k-')
xlabel('$k \: \rm (m)$')
ylabel('$\omega \: \rm (s^{-1})$')
legend('Data','$\omega = \sqrt{gk}$','Location','southeast')
set(findall(gcf,'-property','FontSize'),'FontSize',13)

figname = [fig_folder 'Dipersion_relation'];
saveas(gcf,figname,'fig')

%% Get attenuation coefficient 

disp('Getting attenuation coefficient')
selected_freq = [0.15 0.8]; % selected frequencies between which we proceed to the analysis
x_bound = [1 size(Vx,1)]; % selected boundaries at which we perform 2D FFT
freq_thresh = 0.9;
x_bound_thresh = [round(3*facq_x) round(40*facq_x)];

new_folder_fig = [fig_folder 'Attenuation_fit/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

fig_name = 'Decay_law_video';
save_image = 1;
save_video = 1;
left_bool = 0;

[lambda,dist_fit,freq] = get_attenuation_coef(FFT_t,f,facq_x,selected_freq,x_bound,freq_thresh,x_bound_thresh,left_bool,new_folder_fig,fig_name,save_image,save_video);

% save data to plot attenuation coefficient
freq_min_txt = replace(num2str(selected_freq(1)),'.','p');
freq_max_txt = replace(num2str(selected_freq(2)),'.','p');
f_thresh_txt = replace(num2str(freq_thresh),'.','p');
attenuation_file = [fig_folder 'attenuation_coef_data_fmin' freq_min_txt '_fmax' freq_max_txt '_fthresh' f_thresh_txt];
%attenuation_file = replace(attenuation_file,'.','p');
save(attenuation_file,'lambda','dist_fit','freq','selected_freq','x_bound','freq_thresh','x_bound_thresh')

disp('DONE.')

%% Quick plot of the attenuation law 

%  masking according to dist_fit
mask = (dist_fit < 0.15) & (abs(lambda) > 0.001) ; % keep only points for which distance is smaller than..
f = freq';
fitted_f = f(mask);
fitted_alpha = abs(lambda(mask));

save_boolean = 1;

l1 = fminsearch(@(s)powerfit(fitted_f,fitted_alpha,s),[1,1]);
% l1 = fminsearch(@(s)powerfit(fitted_k,fitted_omega,s),[1,1]);
f_list = linspace(0.01,2,100);
yth = powerfun(f_list,l1); % fitted exponential function

% Take more datas 

% theory = powerfun(freq_secondary,l1);
% Need to take the orthogonal distance in the logarithm scale....
% dist_to_fit = sum((log((theory)) - log((lambda_secondary))).^2)/sum(log(abs(lambda_secondary)).^2);
% 
% mask = dist_to_fit < 0.5;
% freq_plot = freq_secondary(mask);
% lambda_plot = lambda_secondary(mask);

attenuation_fig = figure;
% loglog(freq_secondary,lambda_secondary,'o');
% hold on 
loglog(fitted_f,fitted_alpha,'o','MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerEdgeColor','black');
hold on
plot(f_list,yth,'r--','LineWidth',1.5);
% plot(freq_secondary,theory,'k-')
xlabel('$f \: \rm (Hz)$','Interpreter','latex');
ylabel('$\alpha \: \rm (m^{-1})$','Interpreter','latex');
grid on 
% axis([0.05 2 0.0001 1])
ax = gca;
ax.FontSize = 13;

power_law_txt = ['$\alpha(f) = ' sprintf('%0.2f',l1(2)) 'f^{' sprintf('%0.2f',l1(1)) '}$'];
% legend('Data','Fitted Data',power_law_txt,'Interpreter','latex','Location','southeast','FontSize',13)
legend('Data',power_law_txt,'Interpreter','latex','Location','northwest','FontSize',13)

set_Papermode(gcf);

if save_boolean
    attenuation_filename = [fig_folder 'attenuation_law'];
    saveas(attenuation_fig,attenuation_filename,'fig');
    saveas(attenuation_fig,attenuation_filename,'pdf');
end 

%%

function yth = powerfun(x,l)
    yth = l(2)*x.^l(1);
end 

function d = powerfit(x,y,l)
    yth = powerfun(x,l);
    d = sum((yth-y).^2);
end 
