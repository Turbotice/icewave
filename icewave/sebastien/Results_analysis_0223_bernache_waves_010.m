%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

%base = 'W:/SagWin2024/Data/0223/Drones/bernache/matData/12-waves_010/';
base = 'H:/Rimouski_2024/Data/2024/0223/bernache/matData/12-waves_010/';
% base = 'E:/Rimouski_2024/Data/2024/0219/matData/waves_012/';

filename = 'PIV_processed_i00_N0_Dt4_b1_W32_xROI1_width3388_yROI1_height2159_scaled.mat';
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

facq_x = 1/m.SCALE.fx; % scale in box / meter
facq_t = 1/m.SCALE.ft; % scale in frame / sec
scale_V = m.SCALE.scale_V; % factor scaling for velocity in meter / s
W = m.PIV_param.w ;
font_size = 13;

%% Get Histogram displacement 

% Vx = m.Vx(10:end,:,:) - mean(mean(m.Vx(10:end,:,:),2),1);
filename = [fig_folder 'histogramm_displacements_Vx_time_average'];
average_bool = 1;
get_histogram_displacement(m.Vx/scale_V,W,average_bool,font_size,filename);

%% Get profile of velocity
i_x = 30;
i_y = 60;
left_bool = 1;
Vx = flip(m.Vx,1);
profile_fig = plot_located_profile(Vx,i_x,i_y,facq_x,facq_t,left_bool);
set_Papermode(profile_fig)

disp(mean(abs(Vx(i_x,i_y,:))))

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
Vx = flip(Vx,1); % waves are initially coming from the right boarder of the image
Vy = flip(Vy,1);

disp('Drone motion corrected')

%% Plot main features of the velocity field
disp('Get velocity features')
idx_frame = 20;
caxis_amp = 1.0; % Amplitude of the colorbar scale (in meter / second)

plot_velocity_features(Vx,Vy,facq_x,idx_frame,caxis_amp,fig_folder,1);


%% computes Quantile of the velocity field 
Q = quantile(Vx,[0.1 0.9],'all');

%%
%fig_folder_E = 'E:/Rimouski_2024/Data/2024/0226/Drones/mesange/12-FRAC_001/';
caxis_amp = [Q(1)-0.1 Q(2)+0.1]; 
fig_name = [fig_folder 'Velocity_field_Vx_movie_test'];
fps = facq_t ;
step = 5; % prints every 5 frames 
movie_velocity_field(m.Vx,facq_x,facq_t,caxis_amp,fps,fig_name,0)

%% Show quickly the velocity field 
step = 2; % step between two frames 
fig_name = [fig_folder 'Velocity_field_Vx_movie_step_' num2str(step)];
video_filename = [fig_name '.avi']; % folder where the video is saved
vid = VideoWriter(video_filename);
vid.FrameRate = facq_t/step;
open(vid)

figure, 
x = (1:1:nx)*m.SCALE.fx;
y = (ny:-1:1)*m.SCALE.fx;
t = (1:1:nt)*m.SCALE.ft;
relevant_idx = (1:step:nt);

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
    caxis([Q(1)-0.1 , Q(2)+0.1])
%     set(gca,'XDir','reverse')
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
selected_freq = [0.15 1.0];
x_bound = [1 size(FFT_t,1)];
caxis_amp = 0.1; % amplitude of the colorbar in meter/second
fig_name = ['Demodulated_field_Vx_caxis_' num2str(caxis_amp) 'ms'];

save_image = 1;
save_video = 1;
plot_demodulated_field(FFT_t,f,facq_x,selected_freq,x_bound,caxis_amp,fig_folder,fig_name,save_image,save_video)

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

% #####################################
%% #### E(k,f) plot ###################
% #####################################

% FFT 3D of the velocity field 
N = size(Vx);
add_pow2 = [0, 0, 0]; % additional padding for each dimension 
padding = 2.^(nextpow2(N) + add_pow2); % padding for each dimension

disp('Computing FFT 3D')
FFT = fftn(Vx,padding)/numel(Vx); % FFT 
disp('FFT 3D computed')

kx = 2*pi*facq_x*(-padding(1)/2:padding(1)/2-1)/padding(1);
ky = -2*pi*facq_x*(-padding(2)/2:padding(2)/2-1)/padding(2);

% keep only positive frequencies 
FFT_positive = FFT(:,:,1:padding(3)/2 +1);
FFT_positive(:,:,2:end-1) = 2*FFT_positive(:,:,2:end-1);

f = facq_t*(0:padding(3)/2)/padding(3);

% FFT shift for all dimensions
shift = fftshift(fftshift(FFT_positive,2),1);
disp('FFT shifted')

% Radial average for each frequency

radial_step = 1;
% center of the 2D-FFT
x0 = padding(1)/2 + 1;
y0 = padding(2)/2 + 1;

%loop over all frequencies 
for i0 = 1:size(shift,3)
    % demodulated 2D-FFT
    img = shift(:,:,i0);

    [R_tics, R_profile] = radialavg2(abs(img),radial_step,x0,y0);
    E(i0,:) = R_profile; % Amplitude for value (f,k)
end 
disp('Amplitude computed')

k = 2*pi*R_tics*facq_x/padding(1);

%% Plot A(f,k)
omega = 2*pi*f;

g = 9.81; % gravity intensity 
k_list = linspace(0.1,7,100);
deep_water = sqrt(g*k_list);
h_w = 3.8; % water depth
shallow = sqrt(g*h_w*k_list.^2);
yth = sqrt(g*k_list.*tanh(h_w*k_list));

harmonic2 = sqrt(2*g*k_list.*tanh(h_w*k_list/2));
harmonic3 = sqrt(3*g*k_list.*tanh(h_w*k_list/3));

% plot A(omega,k)
fig_Afk = figure; 
pcolor(k,omega,E)
shading interp
xlabel('$k \: \rm (m^{-1})$')
ylabel('$\omega \: \rm (rad.s^{-1})$')
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'ColorScale','log')
% colormap(slanCM('gnuplot2'))
cbar = colorbar();
cbar.Label.String = '$|\hat{V}_x|(k,\omega)$';
cbar.Label.Interpreter = 'latex';
axis([0.1 7, 2*pi*0.15 2*pi*2.2])
hold on 
loglog(k_list,yth,'w--')
hold on
loglog(k_list,harmonic2,'w--')
hold on 
loglog(k_list,harmonic3,'w--')
set_Papermode(gcf)
set(gca,'FontSize',13)
clim([4e-4 2e-2])
% xticks([1e-1 2e-1 3e-1 1 2 3])
% lgnd = legend('',['$\omega^2 = gk \tanh(gh_w) \: h_w = ' num2str(h_w) '\: \rm m$'],'','');
% set(lgnd,'Location','southeast')
% set(lgnd,'color','none')
% loglog(k_list,harmonic2,'w--')
% hold on 
% loglog(k_list,harmonic3,'w--')
% colormap(slanCM('thermal-2'))

fig_filename = ['A_fk_0223_waves_010_with_harmonics_hw' num2str(h_w)];
fig_filename = replace(fig_filename,'.','p');

saveas(fig_Afk,[fig_folder fig_filename],'fig')
saveas(fig_Afk,[fig_folder fig_filename],'pdf')

%% Save data of plots
data_filename = ['Data_A_fk_0223_waves_010'];
data_filename = replace(data_filename,'.','p');

save([fig_folder data_filename],'f','omega','k','E','shift','-v7.3')
disp('Data saved')

%% Load Data for A(f,k) plot
filename = 'Data_A_fk_0223_waves_010.mat' ;

disp('Loading data..')
load([fig_folder filename])
disp('Data loaded')

%% Detect each branch of bound waves

% Select a profile of the A(f,k) : 
selected_freq = 0.63;
[~,idx_freq] = min(abs(selected_freq - f));

% figure(10),
% plot(2*pi./k,E(idx_freq,:))
% xlabel('$\lambda \: \rm (m)$')
% ylabel('$|\hat{V}_x|(k) \: \rm (m.s^{-1})$')
% xlim([0 50])
freq_end = 0.9; % stop loop at this frequency
[~,idx_end] = min(abs(freq_end - f));

figure(12)
for i = 1 : idx_end

    disp(['Omega = ' num2str(omega(i))])
    max_absolute = max(E(i,:));
    profile = E(i,:)/max_absolute;
%     figure(11),
%     plot(k,profile)
%     xlabel('$k \: \rm (m^{-1})$')
%     ylabel('$|\hat{V}_x|(k) \: \rm (m.s^{-1})$')
%     xlim(2*pi*[1/50 1/2])
% 
    [pks,locs,w,prom] = findpeaks(profile,'MinPeakProminence',1.5e-1);
    [y_max,k_peaks] = subpix_precision(profile,k',locs);
    
    M_peaks(i).k = k_peaks;
    M_peaks(i).A = y_max*max_absolute;
    M_peaks(i).width = w;
    M_peaks(i).omega = omega(i);
    
%     findpeaks(profile,k,'MinPeakProminence',1e-1,'Annotate','extents')
    plot(k,profile,'o')
    hold on 
    plot(M_peaks(i).k,y_max,'rd')
    xlim(2*pi*[1/50 1/2])
    hold off
    xlabel('$k \: \rm (rad.m^{-1})$')
    ylabel('$\hat{V_x}(k,f) \: \rm (m.s^{-1})$')
    pause(0.1)
end 

%% Plot detected peaks
c = [0 0.4470 0.7410];

figure(13)
for i = 1 : size(M_peaks,2)
%     disp(M_peaks(i).k)
%     disp(M_peaks(i).omega .* ones(1,length(M_peaks(i).k)))
    plot(M_peaks(i).k,M_peaks(i).omega .* ones(1,length(M_peaks(i).k)),'o','MarkerFaceColor',c,'MarkerEdgeColor','k')
    hold on 
end

xlabel('$k \: \rm (m^{-1})$')
ylabel('$\omega \: \rm (rad.s^{-1})$')

%% Filtering data & creating two arrays

omega_array = [];
k_array = [];
amplitude_array = [];
for i = 1 : size(M_peaks,2)
    current_omeg = M_peaks(i).omega .* ones(1,length(M_peaks(i).k)); 
    current_k = M_peaks(i).k;
    current_A = M_peaks(i).A;
    omega_array = cat(2,omega_array,current_omeg); % store omega 
    k_array = cat(2,k_array,current_k);
    amplitude_array = cat(2,amplitude_array,current_A);
end 

%% 

mask = ~(((omega_array > 3.0) & (k_array < 0.7)) | (omega_array < 1.0)); 
omega_array = omega_array(mask);
k_array = k_array(mask);
A_array = amplitude_array(mask);


mask2 = ((omega_array > 1.37) & (k_array < 0.2)) | ((k_array > 0.45) & (omega_array < 1.35))...
    | ((omega_array > 4.95) & (k_array < 0.89)) | ((omega_array > 5.44) & (k_array < 1.3));
mask2 = ~mask2;
omega_array = omega_array(mask2);
k_array = k_array(mask2);
A_array = A_array(mask2);

figure, 
plot(k_array,omega_array,'o')
grid on 


%%

g = 9.81; % gravity intensity 
k_list = linspace(0.05,6,100);
deep_water = sqrt(g*k_list);
h_w = 3.8; % water depth
shallow = sqrt(g*h_w*k_list.^2);
yth = sqrt(g*k_list.*tanh(h_w*k_list));

harmonic1 = sqrt(1*g*k_list.*tanh(h_w*k_list));
harmonic2 = sqrt(2*g*k_list.*tanh(h_w*k_list/2));
harmonic3 = sqrt(3*g*k_list.*tanh(h_w*k_list/3));

figure(14)
loglog(k_array,omega_array,'o','MarkerFaceColor',c,'MarkerEdgeColor','k')
xlabel('$k \: \rm (m^{-1})$')
ylabel('$\omega \: \rm (rad.s^{-1})$')
grid on 
hold on 
loglog(k_list,harmonic1,'--b')
hold on 
loglog(k_list,harmonic2,'--r')
% hold on 
% loglog(k_list,harmonic3,'--g')

% axis([0.05 4 , 1e-1 30])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;

% figname = [fig_folder 'Harmonics_shallow_water_hw_' replace(num2str(h_w),'.','p') ];
% saveas(gcf,figname,'fig')
% saveas(gcf,figname,'pdf')
% saveas(gcf,figname,'png')

%% Displace harmonics to main branch

% compute a distance to each harmonics 
N = [1,2,3]';
% k = [1 1.5 2];

% #idx1 : N, idx2 : k
omegaN = bound_harmonicN(k_array,N,h_w);
D = sqrt((repmat(omega_array,[3,1]) - omegaN).*2); 
[~,closest_harmonic] = min(D,[],1);

% set lowest values of k to closest_harmonic = 1
closest_harmonic(k_array < 0.5) = 1;

% set lowest values of harmonic2 to closest_harmonic = 2
% figure,
% mask = (k_array < 0.62) & (closest_harmonic == 3);
% loglog(k_array(mask),omega_array(mask),'o')
% axis([1e-1 2, 5e-1 5])
% closest_harmonic((k_array < 0.62) & (closest_harmonic == 3)) = 2;

%% Plot harmonics with different colours 

% colors = ["blue","red","green"];
% colors = ["#0072BD","#D95319","#77AC30"];
colors = ["#0072BD","#D95319"];
figure(15)
for i = 1 :2
    current_omega = omega_array(closest_harmonic == i);
    current_k = k_array(closest_harmonic == i);
    loglog(current_k,current_omega,'o','MarkerEdgeColor',colors(i))
    hold on
    
end 

xlabel('$k \: \rm (m^{-1})$')
ylabel('$\omega \: \rm (rad.s^{-1})$')
grid on 
hold on 
loglog(k_list,harmonic1,'--b')
hold on 
loglog(k_list,harmonic2,'--r')
% hold on 
% loglog(k_list,harmonic3,'--g')
% axis([0.05 4 , 2e-1 10])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;
% 
% figname = [fig_folder 'Harmonics_shallow_water_hw_' replace(num2str(h_w),'.','p') '_different_colors'];
% saveas(gcf,figname,'fig')
% saveas(gcf,figname,'pdf')
% saveas(gcf,figname,'png')

%% Save selected points

filename = 'Data_plot_selected_harmonics_0223_bernache_waves_010';
directory = 'E:/Rimouski_2024/Data/2024/0223/bernache/matData/12-waves_010/Plots/';
save([directory filename],'omega_array','k_array','A_array','closest_harmonic','colors','k_list','harmonic1','harmonic2','-v7.3')

%% Plot each harmonic with a different colormaps corresponding to the intensity of each points (omega,k)

% Load Data of selected harmonics 
filename = 'Data_plot_selected_harmonics_0223_bernache_waves_010';
base = 'H:/Rimouski_2024/Data/2024/0223/bernache/matData/12-waves_010/Plots/';
load([base filename])
disp('Data of separated harmonics loaded')

%%
A_h1 = A_array(closest_harmonic == 1);
A_h1 = (A_h1 - min(A_h1))/(max(A_h1) - min(A_h1));

A_h2 = A_array(closest_harmonic == 2);
A_h2 = (A_h2 - min(A_h2))/(max(A_h2) - min(A_h2));
marker_size = 50;

harmonic_cmaps = figure;
% create two different axes 
ax1 = axes(harmonic_cmaps);
ax2 = copyobj(ax1,harmonic_cmaps);

plot(ax1,k_list,harmonic1,'--b')
hold on
plot(ax1,k_list,harmonic2,'--r')
hold on 
s1 = scatter(ax1,k_array(closest_harmonic == 1),omega_array(closest_harmonic == 1),marker_size,A_h1,'filled','MarkerEdgeColor','k');

s2 = scatter(ax2,k_array(closest_harmonic == 2),omega_array(closest_harmonic == 2),marker_size,A_h2,'filled','MarkerEdgeColor','k');

colormap(ax1,slanCM(9))
colormap(ax2,slanCM(11))
ax2.Visible = "off";

% sets limits of each axes
xlim(ax1,[0.1 5])
xlim(ax2,[0.1 5])
ylim(ax1,[0.7 8])
ylim(ax2,[0.7 8])
grid on 
% set axes scale in log-log
set(ax1,'xscale','log')
set(ax1,'yscale','log')
set(ax2,'xscale','log')
set(ax2,'yscale','log')

xlabel(ax1,'$k \: \rm (rad.m^{-1})$')
ylabel(ax1,'$\omega \: \rm (rad.s^{-1})$')
legend(ax1,['$\omega_1 = \sqrt{gk \tanh(' num2str(h_w) 'k)}$'],...
    ['$\omega_2 = \sqrt{2gk \tanh(' num2str(h_w) '\frac{k}{2})}$'],'','Location','southeast')
ax1.FontSize = 13;
set_Papermode(harmonic_cmaps)

fig_filename = [fig_folder 'Harmonics_colormaps_12-waves_010_hw_' replace(num2str(h_w),'.','p') ];
saveas(harmonic_cmaps,fig_filename,'fig')
saveas(harmonic_cmaps,fig_filename,'pdf')

%% Move each branch back to main one 

harmonic_cmaps = figure;
% create two different axes 
ax1 = axes(harmonic_cmaps);
ax2 = copyobj(ax1,harmonic_cmaps);

plot(ax1,k_list,harmonic1,'--b')
hold on
plot(ax1,k_list,harmonic2,'--r')
hold on 
s1 = scatter(ax1,k_array(closest_harmonic == 1),omega_array(closest_harmonic == 1),marker_size,A_h1,'filled','MarkerEdgeColor','k');

s2 = scatter(ax2,k_array(closest_harmonic == 2)/2,omega_array(closest_harmonic == 2)/2,marker_size,A_h2,'filled','MarkerEdgeColor','k');

colormap(ax1,slanCM(9))
colormap(ax2,slanCM(11))
ax2.Visible = "off";

% sets limits of each axes
xlim(ax1,[0.1 5])
xlim(ax2,[0.1 5])
ylim(ax1,[0.7 8])
ylim(ax2,[0.7 8])
grid on 
% set axes scale in log-log
set(ax1,'xscale','log')
set(ax1,'yscale','log')
set(ax2,'xscale','log')
set(ax2,'yscale','log')

xlabel(ax1,'$k \: \rm (rad.m^{-1})$')
ylabel(ax1,'$\omega \: \rm (rad.s^{-1})$')
legend(ax1,['$\omega_1 = \sqrt{gk \tanh(' num2str(h_w) 'k)}$'],...
    ['$\omega_2 = \sqrt{2gk \tanh(' num2str(h_w) '\frac{k}{2})}$'],'','Location','southeast')

set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;
figname = [fig_folder 'Recomposition_harmonics_water_12-waves_010_hw_' replace(num2str(h_w),'.','p')];
saveas(gcf,figname,'fig')
saveas(gcf,figname,'pdf')
saveas(gcf,figname,'png')

%%
figure(16)
for i = 1 :2
    
    current_omega = omega_array(closest_harmonic == i);
    current_k = k_array(closest_harmonic == i);
    if i > 1
        current_omega = current_omega/i;
        current_k = current_k/i;
    end 
    loglog(current_k,current_omega,'o','MarkerFaceColor',colors(i),'MarkerEdgeColor',colors(i),'MarkerSize',5)
    hold on
    
end 

xlabel('$k \: \rm (m^{-1})$')
ylabel('$\omega \: \rm (rad.s^{-1})$')
grid on 
hold on 
loglog(k_list,harmonic1,'--b')
hold on 
loglog(k_list,harmonic2,'--r')
% hold on 
% loglog(k_list,harmonic3,'--g')
axis([0.1 6 , 5e-1 10])


% #######################################
%% Attenuation coefficient in corrected direction 
% #######################################

% Parameters 
selected_freq = [0.3 0.6]; % selected frequencies between which we proceed to the analysis
% get indices of frequencies closest to max and min selected frequency
[min_freq, start_freq_idx] = min(abs(f - selected_freq(1))); 
[max_freq, end_freq_idx] = min(abs(f - selected_freq(2)));
new_freq = f(start_freq_idx : end_freq_idx); % new frequency array 
FFT_cropped = FFT_t(:,:,start_freq_idx:end_freq_idx);
nf = length(new_freq);

% Select range of frequencies and space of shortened fit 
f_short = [0.45 0.54]; % frequency above which the attenuation fit is shortened
x_short = [50 30]; % meters % range over which the shrotened attenuation fit is performed

x = m.x;
y = m.y;
padding_bool = 1;
add_pow2 = 2;
threshold = 0.8;
L0 = 80; % segments size in meter 

new_folder_fig = [fig_folder 'attenuation_oriented/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

alpha = zeros(nf,1); % array of attenuation coefficients
C  = zeros(nx,1);
d = zeros(nf,1); % array of distance to plot

for idx_freq = 1:nf
    current_freq = new_freq(idx_freq);
    disp(['f = ' num2str(current_freq) ' Hz'])
    field = FFT_cropped(:,:,idx_freq);
    freq_txt = replace(num2str(new_freq(idx_freq)),'.','p');
    % 2D FFT of this complex field 
    [shifted_fft,fft_2D,kx,ky] = spatial_FFT(field,padding_bool,add_pow2,facq_x);

    % Detect peaks 
    zmax = max(abs(shifted_fft),[],'all');
    norm_fft = abs(shifted_fft)/zmax;

    peak_idx = 1;
    % binarize the 2D-FFT spectrum
    bin_fft = norm_fft;
    bin_fft(norm_fft > threshold) = 1;
    bin_fft(norm_fft <= threshold) = 0;

    CC = bwconncomp(bin_fft');
    stats = regionprops("table",CC,norm_fft','WeightedCentroid');
    row = round(stats.WeightedCentroid(peak_idx,1));
    col = round(stats.WeightedCentroid(peak_idx,2));
    kx_peak = kx(row);
    ky_peak = ky(col);

    FFT_space_fig = figure(1);
    pcolor(kx,ky,abs(shifted_fft)')
    shading interp
    hold on 
    plot(kx_peak,ky_peak,'ro')
    xlabel('$k_x \: \rm (rad.m^{-1})$')
    ylabel('$k_y \: \rm (rad.m^{-1})$')
    axis([-2 2 -2 2 ]) 
    hold off

    % select peak index 
    % prompt = "Which peak index would you choose? (sorted in amplitude) : ";
    % peak_idx = str2double(input(prompt,"s"));
    % 
    % row = round(stats.WeightedCentroid(peak_idx,1));
    % col = round(stats.WeightedCentroid(peak_idx,2));
    % kx_peak = kx(row);
    % ky_peak = ky(col);
    
    [theta,k_peak] = cart2pol(kx_peak,ky_peak); % converts to polar coordinates 
    disp(['Theta = ' num2str(theta)])

    if theta < 0
        % Build several segments for theta < 0
        % initial line 
        x0 = m.x(1); % minimal value along x axis
        y0 = L0*sin(abs(theta)) + m.y(1);
        ds = m.SCALE.fx; % step of curvilinear coordinate
        s = (0:ds:L0); % curvilinear coordinate
        
        % Points that define the line 
        x_line = x0+s*cos(theta);
        y_line = y0+s*sin(theta);

    else
        % Building several segments for theta > 0
        x0 = (m.y(end) - L0*sin(theta) - m.y(1))*tan(theta);
        y0 = m.y(1);
        ds = m.SCALE.fx;

        s = (0:ds:L0); % curvilinear coordinate
        x_line = x0 + s*cos(theta);
        y_line = y0 + s*sin(theta);
    end 

    real_field_fig = figure(2); 
    pcolor(x,y,real(field)')
    shading interp
    colormap(redblue)
    hold on 
    plot(x_line,y_line,'k--')
    axis image
    xlabel('$x \: \rm (m)$')
    ylabel('$y \: \rm (m)$')
    ax = gca; 
    ax.FontSize = 13;
    set_Papermode(gcf)
 
    % ask if we choose reoriented field
    prompt = "Do you want to reorient wave field ? y/n [y]: ";
    txt = input(prompt,"s");
    if isempty(txt)
        txt = 'y';
    end

    if txt == 'y'
        % Compute field correctly oriented
        [field_star,x_star,y_star] = orientation_parallel2propagation(field,theta,m.SCALE.fx,L0);
        hold on 
        plot(x_line,y_line,'k-')
    else
        field_star = field;
        x_star = m.x;
        y_star = m.y;
    end
      
    hold off
    figname = [new_folder_fig 'Real_field_f' freq_txt];
    saveas(gcf,figname,'fig')
    saveas(gcf,figname,'pdf')

    % Exponential fit 
    TF_x = squeeze(mean(field_star,2)); % average the TF over y  
    A = abs(TF_x); % amplitude along x-axis for each frequency
    
    % indices used to fit the exponential decay

    if isempty(f_short)
        i_min = 1;
        i_max = length(A);
    else 
        for idx_subdomain = 1:length(f_short)
            current_fshort = f_short(idx_subdomain);
            current_xshort = x_short(idx_subdomain);
            [~,i_short] = min(abs(m.x - current_xshort));
            if current_freq < current_fshort 
                i_min = 1;
                i_max = length(A);
            else
                i_min = 1;
                i_max = i_short;
            end 
        end 
    end 

    decay_fig = figure(3);
    decay_fig.Color = [1,1,1];

    log_A = log10(A); % take the log10 of the amplitude of freq i
    % restrict to the boundaries in order to fit
    x_fit = x_star(i_min:i_max); % in meter !!
    A_fit = log_A(i_min:i_max); 
    p = polyfit(x_fit,A_fit,1); % fit log_A by a degree 1 polynome
    alpha(idx_freq) = log(10)*p(1); % get attenuation coefficient in m^-1
    C(idx_freq) = 10.^p(2); % prefactor

    disp(['alpha = ' num2str(alpha(idx_freq)) ' m-1'])
    y_poly = 10.^polyval(p,x_fit);
    plot(x_fit,A(i_min:i_max),'o');
    hold on 
    plot(x_fit,y_poly,'r');
    xlabel('$x \: \rm (m)$','Interpreter','latex');
    ylabel('$\langle | \hat{V_x} | \rangle _y (x,f) \: \rm (m)$','Interpreter','latex');
    grid on 
    data_txt = 'Data';
    fit_txt = ['$y(x) = ' sprintf('%0.2f',C(idx_freq)) ' e^{' sprintf('%0.3f',alpha(idx_freq)) 'x}$'];
    legend(data_txt,fit_txt,'Interpreter','latex','location','northeast','FontSize',13);
    ax = gca;
    ax.FontSize = 13;
    set_Papermode(decay_fig);
    hold off
    
    figname = [new_folder_fig 'Attenuation' freq_txt];
    saveas(decay_fig,figname,'fig')
    saveas(decay_fig,figname,'pdf')
    
    A_red = A(i_min:i_max); % restricted to the region of interest
    d(idx_freq) = sum((y_poly - A_red').^2)/sum(A_red.^2); % distance to the fit

end 

%% Save relevant variables 
attenuation_file =  ['Data_attenuation_oriented_' num2str(selected_freq(1)) 'Hz_to_' num2str(selected_freq(2)) 'Hz'];
attenuation_file = replace(attenuation_file,'.','p');
attenuation_file = [fig_folder attenuation_file '_mode_1_analyse_1021'];
save(attenuation_file,'alpha','C','d','new_freq','selected_freq','L0','threshold','padding_bool','add_pow2','f_short','x_short')

disp('DONE.')

% ##################################
%% Attenuation coefficient mode 1 
% ##################################

% load parameters 
filename = attenuation_file;
S = load(filename);
disp('Attenuation data loaded')

%%

%  masking according to dist_fit
d_thresh = 0.30;

mask = (S.d < d_thresh) & (abs(S.alpha) > 0.001) ; % keep only points for which distance is smaller than..
f = S.new_freq;
fitted_f = f(mask);
fitted_alpha = abs(S.alpha(mask))';

l1 = fminsearch(@(s)powerfit(fitted_f,fitted_alpha,s),[1,1]);
f_list = linspace(0.01,10,100);
yth = powerfun(f_list,l1); % fitted exponential function

attenuation_fig = figure;
loglog(fitted_f,fitted_alpha,'o','MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerEdgeColor','black');
hold on
plot(f_list,yth,'r--','LineWidth',1.5);
xlabel('$f \: \rm (Hz)$','Interpreter','latex');
ylabel('$\alpha \: \rm (m^{-1})$','Interpreter','latex');
grid on 
axis([4e-2 4 1e-3 1])
ax = gca;
ax.FontSize = 13;

power_law_txt = ['$\alpha(f) = ' sprintf('%0.2f',l1(2)) 'f^{' sprintf('%0.2f',l1(1)) '}$'];
% legend('Data','Fitted Data',power_law_txt,'Interpreter','latex','Location','southeast','FontSize',13)
legend('',power_law_txt,'Interpreter','latex','Location','northwest','FontSize',13)
set_Papermode(gcf);

thresh_txt = replace(num2str(1-d_thresh),'.','p');
attenuation_filename = [fig_folder 'attenuation_law_mode_1_confidence_' thresh_txt];
saveas(attenuation_fig,attenuation_filename,'fig');
saveas(attenuation_fig,attenuation_filename,'pdf');

%% Plot distance d to the curve as function of frequencies
figure, 

plot(S.new_freq,S.d,'o')
xlabel('$f \: \rm (Hz)$')
ylabel('$d$')

% ####################
%% MAIN DEVELOPMENTS 
% ####################

% Check reorientation of wave field 
% Parameters 
selected_freq = 0.42; % selected frequencies between which we proceed to the analysis
[~,idx_freq] = min(abs(f- selected_freq));

field  = FFT_t(:,:,idx_freq);
padding_bool = 1;
add_pow2 = 2;
threshold = 0.8;
L0 = 80; % segments size in meter 
[shifted_fft,fft_2D,kx,ky] = spatial_FFT(field,padding_bool,add_pow2,facq_x);

% Detect peaks 
zmax = max(abs(shifted_fft),[],'all');
norm_fft = abs(shifted_fft)/zmax;

i = 1; % select first peak 
% binarize the 2D-FFT spectrum
bin_fft = norm_fft;
bin_fft(norm_fft > threshold) = 1;
bin_fft(norm_fft <= threshold) = 0;

CC = bwconncomp(bin_fft');
stats = regionprops("table",CC,norm_fft','WeightedCentroid');
row = round(stats.WeightedCentroid(i,1));
col = round(stats.WeightedCentroid(i,2));
kx_peak = kx(row);
ky_peak = ky(col);

[theta,k_peak] = cart2pol(kx_peak,ky_peak); % converts to polar coordinates 
disp(['Theta = ' num2str(theta)])

if theta < 0
    % Build several segments for theta < 0
    % initial line 
    x0 = m.x(1); % minimal value along x axis
    y0 = L0*sin(abs(theta)) + m.y(1);
    ds = m.SCALE.fx; % step of curvilinear coordinate
    s = (0:ds:L0); % curvilinear coordinate
    
    % Points that define the line 
    x_line = x0+s*cos(theta);
    y_line = y0+s*sin(theta);

    Nb_lines = floor((y(end) - y0)/(ds*cos(theta))); % number of lines to draw 
    X_line = zeros(length(s),Nb_lines);
    Y_line = zeros(length(s),Nb_lines);

    for j = 1:Nb_lines
        x0 = x0 + ds*sin(abs(theta));
        y0 = y0 + ds*cos(theta);

        X_line(:,j) = x0 + s*cos(theta); % #1 : s-coordinate, #2 line index 
        Y_line(:,j) = y0 + s*sin(theta);
    end
else 
    x0 = (m.y(end) - L0*sin(theta) - m.y(1))*tan(theta);
    y0 = m.y(1);
    ds = m.SCALE.fx;

    s = (0:ds:L0); % curvilinear coordinate
    x_line = x0 + s*cos(theta);
    y_line = y0 + s*sin(theta);

    Nb_lines = floor((y(end) - L0*sin(theta) - y(1))/(ds*cos(theta))); % number of lines to draw 
    X_line = zeros(length(s),Nb_lines);
    Y_line = zeros(length(s),Nb_lines);
    
    for j = 1:Nb_lines
        x0 = x0 - ds*sin(abs(theta));
        y0 = y0 + ds*cos(theta);

        X_line(:,j) = x0 + s*cos(theta); % #1 : s-coordinate, #2 line index 
        Y_line(:,j) = y0 + s*sin(theta);
    end
end 


figure(1), 
pcolor(m.x,m.y,real(field)')
shading interp
colormap(redblue)
hold on 
for j = 1 : size(X_line,2)
    plot(X_line(:,j),Y_line(:,j),'k--')
end 
axis image
xlabel('$x \: \rm (m)$')
ylabel('$y \: \rm (m)$')
ax = gca; 
ax.FontSize = 13;

figure(2),
pcolor(kx,ky,abs(shifted_fft)')
shading interp
axis image
xlabel('$k_x \: \rm (rad.m^{-1})$')
ylabel('$k_y \: \rm (rad.m^{-1})$')


%%
    
    real_field_fig = figure(1); 
    pcolor(x,y,real(field)')
    shading interp
    colormap(redblue)
    hold on 
    plot(x_line,y_line,'k--')
    axis image
    xlabel('$x \: \rm (m)$')
    ylabel('$y \: \rm (m)$')
    ax = gca; 
    ax.FontSize = 13;
    set_Papermode(gcf)
    hold off
    
    figname = [new_folder_fig 'Real_field_f' freq_txt];
    saveas(gcf,figname,'fig')
    saveas(gcf,figname,'pdf')
    
    pause(0.2)
    % Compute field correctly oriented
    [field_star,x_star,y_star] = orientation_parallel2propagation(field,theta,m.SCALE.fx,L0);
    


%% Scale an image of the video 
% Initial picture with scaling 
picture_file = 'W:/SagWin2024/Data/0226/Drones/mesange/matData/10-waves_005/im_0000.tiff';
img = imread(picture_file);

%%
% x = (1:1:size(img,2));
% y = (1:1:size(img,1));

font_size = 13;
new_img = img(:,600:3240,:);
RI = imref2d(size(new_img));

RI.XWorldLimits = RI.XWorldLimits ./m.facq_pix; 
RI.YWorldLimits = RI.YWorldLimits ./m.facq_pix;
img_fig = figure; 
imshow(new_img,RI)
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
axis image
ax = gca;
ax.FontSize = font_size;
axis('xy')
% set correctly the image position for a pdf format 
set_Papermode(gcf)

figname = [fig_folder 'Im_0000_ice'];
saveas(img_fig,figname,'fig')
saveas(img_fig,figname,'pdf')

%% Previous versions 


alpha = zeros(nf,1); % array of attenuation coefficients
d = zeros(nf,1); % array of distance to plot

for idx_freq = 1:nf
    current_freq = new_freq(idx_freq);
    disp(['f = ' num2str(current_freq) ' Hz'])
    field = FFT_cropped(:,:,idx_freq);
    freq_txt = replace(num2str(new_freq(idx_freq)),'.','p');
    % 2D FFT of this complex field 
    [shifted_fft,fft_2D,kx,ky] = spatial_FFT(field,padding_bool,add_pow2,facq_x);

    % Detect peaks 
    zmax = max(abs(shifted_fft),[],'all');
    norm_fft = abs(shifted_fft)/zmax;

    peak_idx = 1;
    % binarize the 2D-FFT spectrum
    bin_fft = norm_fft;
    bin_fft(norm_fft > threshold) = 1;
    bin_fft(norm_fft <= threshold) = 0;

    CC = bwconncomp(bin_fft');
    stats = regionprops("table",CC,norm_fft','WeightedCentroid');
    row = round(stats.WeightedCentroid(peak_idx,1));
    col = round(stats.WeightedCentroid(peak_idx,2));
    kx_peak = kx(row);
    ky_peak = ky(col);

    FFT_space_fig = figure(1);
    pcolor(kx,ky,abs(shifted_fft)')
    shading interp
    hold on 
    plot(kx_peak,ky_peak,'ro')
    xlabel('$k_x \: \rm (rad.m^{-1})$')
    ylabel('$k_y \: \rm (rad.m^{-1})$')
    axis([-2 2 -2 2 ]) 
    hold off

    % select peak index 
    % prompt = "Which peak index would you choose? (sorted in amplitude) : ";
    % peak_idx = str2double(input(prompt,"s"));
    % 
    % row = round(stats.WeightedCentroid(peak_idx,1));
    % col = round(stats.WeightedCentroid(peak_idx,2));
    % kx_peak = kx(row);
    % ky_peak = ky(col);
    % 
    % [theta,k_peak] = cart2pol(kx_peak,ky_peak); % converts to polar coordinates 
    % disp(['Theta = ' num2str(theta)])

    if theta < 0
        % Build several segments for theta < 0
        % initial line 
        x0 = m.x(1); % minimal value along x axis
        y0 = L0*sin(abs(theta)) + m.y(1);
        ds = m.SCALE.fx; % step of curvilinear coordinate
        s = (0:ds:L0); % curvilinear coordinate
        
        % Points that define the line 
        x_line = x0+s*cos(theta);
        y_line = y0+s*sin(theta);

        % Compute field correctly oriented
        [field_star,x_star,y_star] = orientation_parallel2propagation(field,theta,m.SCALE.fx,L0);

    else
        % Building several segments for theta > 0
        x0 = (m.y(end) - L0*sin(theta) - m.y(1))*tan(theta);
        y0 = m.y(1);
        ds = m.SCALE.fx;

        s = (0:ds:L0); % curvilinear coordinate
        x_line = x0 + s*cos(theta);
        y_line = y0 + s*sin(theta);
    end 

    real_field_fig = figure(2); 
    pcolor(x,y,real(field)')
    shading interp
    colormap(redblue)
    hold on 
    plot(x_line,y_line,'k--')
    axis image
    xlabel('$x \: \rm (m)$')
    ylabel('$y \: \rm (m)$')
    ax = gca; 
    ax.FontSize = 13;
    set_Papermode(gcf)
 
    % ask if we choose reoriented field
    prompt = "Do you want to reorient wave field ? y/n [y]: ";
    txt = input(prompt,"s");
    if isempty(txt)
        txt = 'y';
    end

    if txt == 'y'
        % Compute field correctly oriented
        [field_star,x_star,y_star] = orientation_parallel2propagation(field,theta,m.SCALE.fx,L0);
        hold on 
        plot(x_line,y_line,'k-')
    else
        field_star = field;
        x_star = m.x;
        y_star = m.y;
    end
      
    hold off
    figname = [new_folder_fig 'Real_field_f' freq_txt];
    saveas(gcf,figname,'fig')
    saveas(gcf,figname,'pdf')

    % Exponential fit 
    TF_x = squeeze(mean(field_star,2)); % average the TF over y  
    A = abs(TF_x); % amplitude along x-axis for each frequency
    
    % indices used to fit the exponential decay
    for idx_subdomain = 1:length(f_short)
        current_fshort = f_short(idx_subdomain);
        current_xshort = x_short(idx_subdomain);
        [~,i_short] = min(abs(m.x - current_xshort));
        if current_freq < current_fshort 
            i_min = 1;
            i_max = length(A);
        else
            i_min = 1;
            i_max = i_short;
        end 
    end 

    decay_fig = figure(3);
    decay_fig.Color = [1,1,1];

    log_A = log10(A); % take the log10 of the amplitude of freq i
    % restrict to the boundaries in order to fit
    x_fit = x_star(i_min:i_max); % in meter !!
    A_fit = log_A(i_min:i_max); 
    p = polyfit(x_fit,A_fit,1); % fit log_A by a degree 1 polynome
    alpha(idx_freq) = log(10)*p(1); % get attenuation coefficient in m^-1

    disp(['alpha = ' num2str(alpha(idx_freq)) ' m-1'])
    y_poly = 10.^polyval(p,x_fit);
    plot(x_fit,A(i_min:i_max),'o');
    hold on 
    plot(x_fit,y_poly,'r');
    xlabel('$x \: \rm (m)$','Interpreter','latex');
    ylabel('$\langle | \hat{V_x} | \rangle _y (x,f) \: \rm (m)$','Interpreter','latex');
    grid on 
    data_txt = 'Data';
    fit_txt = ['$y(x) = C e^{' sprintf('%0.3f',alpha(idx_freq)) 'x}$'];
    legend(data_txt,fit_txt,'Interpreter','latex','location','northeast','FontSize',13);
    ax = gca;
    ax.FontSize = 13;
    set_Papermode(decay_fig);
    hold off
    
    figname = [new_folder_fig 'Attenuation' freq_txt];
    saveas(decay_fig,figname,'fig')
    saveas(decay_fig,figname,'pdf')
    
    A_red = A(i_min:i_max); % restricted to the region of interest
    d(idx_freq) = sum((y_poly - A_red').^2)/sum(A_red.^2); % distance to the fit

end 









%%

function yth = powerfun(x,l)
    yth = l(2)*x.^l(1);
end 

function d = powerfit(x,y,l)
    yth = powerfun(x,l);
    d = sum((yth-y).^2);
end 


function [omegaN] = bound_harmonicN(k,N,h_w)
    % Computes the pulsation omegaN for a given k, along the dispersion
    % relation of the harmonic N
    
    omegaN = sqrt(9.81*N.*k.*tanh(h_w.*k./N));
end 

function yth = lorentzian_fun_2param(x,P)
    % This function enables to fit a curve with a lorentzian function,
    % using only two parameters : P(1) and P(2)
    yth = P(1) ./ (P(2) + x.^2);
end 

function d = lorentzian_fit_2param(x,y,P)

    yth = lorentzian_fun_2param(x,P);
    d = sum((yth - y).^2);

end 
% function [y_max,x_max] = subpix_precision(y,x,j0)
%     % Computes subpixel precision using 2nd Order fitting
%     % Arguments : - y : amplitude array 
%     %             - x : abscisse array 
%     %             - j0 : closest indices to local maxima
% 
%     % Returns : - y_max : maximum of the fitted curve (unit of y)
%     %           - x_max : subpixel value of x for which y_max is reached (unit
%     %           of x)
%     p_0 = y(j0);
%     p_right = y(j0 + 1);
%     p_left = y(j0 - 1);
% 
%     delta = (p_right - p_left)./(2*(2*p_0 - p_right - p_left));
%     x_max = x(j0) + (x(j0+1) - x(j0)).*delta;
%     % Computation of the subpixel amplitude
%     y_max = p_0 + 0.5*(2*p_0 - p_right - p_left).*delta.^2;
% 
% end
