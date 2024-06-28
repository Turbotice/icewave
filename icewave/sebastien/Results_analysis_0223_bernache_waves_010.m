%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

base = 'W:/SagWin2024/Data/0223/Drones/bernache/matData/12-waves_010/';
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

%% Get Histogram displacement 

% Vx = m.Vx(10:end,:,:) - mean(mean(m.Vx(10:end,:,:),2),1);
filename = [fig_folder 'histogramm_displacements_Vx_time_average'];
average_bool = 1;
get_histogram_displacement(m.Vx/scale_V,W,average_bool,font_size,filename);

%% Get profile of velocity
i_x = 30;
i_y = 60;
left_bool = 1;
profile_fig = plot_located_profile(m.Vx,i_x,i_y,facq_x,facq_t,left_bool);
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
disp('Get velocity features')
idx_frame = 20;
caxis_amp = 1.0; % Amplitude of the colorbar scale (in meter / second)

plot_velocity_features(Vx,Vy,facq_x,idx_frame,caxis_amp,fig_folder,1);


%% computes Quantile of the velocity field 
Q = quantile(m.Vx,[0.1 0.9],'all');

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
x = (1:1:nx)*m.fx;
y = (ny:-1:1)*m.fx;
t = (1:1:nt)*m.ft;
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
selected_freq = [0.1 0.8];
x_bound = [1 size(FFT_t,1)];
caxis_amp = 0.02; % amplitude of the colorbar in meter/second
fig_name = ['Demodulated_field_Vx_caxis_' num2str(caxis_amp) 'ms'];

save_image = 1;
save_video = 1;
plot_demodulated_field(FFT_t,f,facq_x,selected_freq,x_bound,caxis_amp,1,fig_folder,fig_name,save_image,save_video)

%% Get wave vectors 
disp('Getting wave vectors')
selected_freq = [0.1 0.8]; % selected frequencies between which we proceed to the analysis
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
% load data 
base_disp = 'W:/SagWin2024/Data/0226/Drones/mesange/matData/10-waves_005/Plots/';
filename_disp = 'dispersion_relation_data_Vx_fmin0p1_fmax0p8_add_pow2.mat';
Sdisp = load([base_disp filename_disp]);
disp('Dispersion relation data loaded')

%%
freq = Sdisp.freq;
omega = 2*pi*freq;
k = Sdisp.k;
mask = (omega > 2.8);
filt_k = k(mask);
filt_omega = omega(mask);

g = 9.81;
k_list = linspace(0.07,5,100);
deep_water = sqrt(g*k_list);

figure,
loglog(k,omega,'o')
hold on 
loglog(filt_k*2,filt_omega,'d')
grid on 
hold on 
loglog(k_list,deep_water,'k-')
xlabel('$k \: \rm (m)$')
ylabel('$\omega \: \rm (rad.s^{-1})$')
legend('Data','sub-k harmonic','$\omega = \sqrt{gk}$','Location','southeast')
set(findall(gcf,'-property','FontSize'),'FontSize',13)

figname = [fig_folder 'Dipersion_relation'];
% saveas(gcf,figname,'fig')

% ##################################
%% Get attenuation coefficient 
% ##################################

disp('Getting attenuation coefficient')
selected_freq = [0.1 0.8]; % selected frequencies between which we proceed to the analysis
x_bound = [1 size(Vx,1)]; % selected boundaries at which we perform 2D FFT
freq_thresh = [0.36  0.39  0.45  0.50  0.55  0.7];
x_bound_thresh = round([1 90; 1 60 ; 1 50; 1 35; 1 20; 1 15]*facq_x);

new_folder_fig = [fig_folder 'Attenuation_fit/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

fig_name = 'Decay_law_video_mode_1';
save_image = 1;
save_video = 1;
left_bool = 1;

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

%% Quick plot of the attenuation law 

figure,
loglog(freq,abs(lambda),'o')
grid on 

%%

%  masking according to dist_fit
d_thresh = 0.05;
mask = (dist_fit < d_thresh) & (abs(lambda) > 0.001) ; % keep only points for which distance is smaller than..
f = freq;
fitted_f = f(mask);
fitted_alpha = abs(lambda(mask))';

save_boolean = 1;

l1 = fminsearch(@(s)powerfit(fitted_f,fitted_alpha,s),[1,1]);
% l1 = fminsearch(@(s)powerfit(fitted_k,fitted_omega,s),[1,1]);
f_list = linspace(0.01,10,100);
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
axis([4e-2 4 1e-3 1])
ax = gca;
ax.FontSize = 13;

power_law_txt = ['$\alpha(f) = ' sprintf('%0.2f',l1(2)) 'f^{' sprintf('%0.2f',l1(1)) '}$'];
% legend('Data','Fitted Data',power_law_txt,'Interpreter','latex','Location','southeast','FontSize',13)
legend('',power_law_txt,'Interpreter','latex','Location','northwest','FontSize',13)

set_Papermode(gcf);
thresh_txt = replace(num2str(1-d_thresh),'.','p');
if save_boolean
    attenuation_filename = [fig_folder 'attenuation_law_mode_1_confidence_' thresh_txt];
    saveas(attenuation_fig,attenuation_filename,'fig');
    saveas(attenuation_fig,attenuation_filename,'pdf');
end 

% #####################################
%% #### E(k,f) plot ###################
% #####################################

% FFT 3D of the velocity field 
N = size(Vx);
add_pow2 = [0, 0, 0]; % additional padding for each dimension 
padding = 2.^(nextpow2(N) + add_pow2); % padding for each dimension

disp('Computing FFT')
FFT = fftn(Vx,padding)/numel(Vx); % FFT 
disp('FFT computed')

%%
kx = 2*pi*facq_x*(-padding(1)/2:padding(1)/2-1)/padding(1);
ky = -2*pi*facq_x*(-padding(2)/2:padding(2)/2-1)/padding(2);

%% keep only positive frequencies 
FFT_positive = FFT(:,:,1:padding(3)/2 +1);
FFT_positive(:,:,2:end-1) = 2*FFT_positive(:,:,2:end-1);

f = facq_t*(0:padding(3)/2)/padding(3);

%% FFT shift for all dimensions
shift = fftshift(fftshift(FFT_positive,2),1);
disp('FFT shifted')

%% Radial average for each frequency
% selected_freq = 0.23; 
% [~,i0] = min(abs(f - selected_freq));

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
k_list = linspace(0.1,6,100);
deep_water = sqrt(g*k_list);
h_w = 3.8; % water depth
shallow = sqrt(g*h_w*k_list.^2);
yth = sqrt(g*k_list.*tanh(h_w*k_list));

harmonic2 = sqrt(2*g*k_list.*tanh(h_w*k_list/2));
harmonic3 = sqrt(3*g*k_list.*tanh(h_w*k_list/3));

% plot A(omega,k)
figure, 
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
axis([0.1 6, 2*pi*0.1 2*pi*2])
% hold on 
% loglog(k_list,yth,'w--')
% hold on
% loglog(k_list,harmonic2,'w--')
% hold on 
% loglog(k_list,harmonic3,'w--')
set_Papermode(gcf)
set(gca,'FontSize',13)
caxis([1e-5 5e-3])
% xticks([1e-1 2e-1 3e-1 1 2 3])
% lgnd = legend('',['$\omega^2 = gk \tanh(gh_w) \: h_w = ' num2str(h_w) '\: \rm m$'],'','');
% set(lgnd,'Location','southeast')
% set(lgnd,'color','none')
% loglog(k_list,harmonic2,'w--')
% hold on 
% loglog(k_list,harmonic3,'w--')
% colormap(slanCM('thermal-2'))

%% Save data of plots
data_filename = ['Data_A_fk_0223_waves_005_hw' num2str(h_w)];
data_filename = replace(data_filename,'.','p');

save([fig_folder data_filename],'f','omega','k','E','shift','-v7.3')
disp('Data saved')

%% Load Data for A(f,k) plot
filename = 'Data_A_fk_0223_waves_005_hw3p8.mat' ;
base = 'W:/SagWin2024/Data/0226/Drones/mesange/matData/10-waves_005/Plots/';

disp('Loading data..')
load([base filename])
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
    [pks,locs,w,prom] = findpeaks(profile,'MinPeakProminence',1e-1);
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

mask = ~(((omega_array > 1.58) & (k_array < 0.28)) | (omega_array < 0.4)) ; 
omega_array = omega_array(mask);
k_array = k_array(mask);
A_array = amplitude_array(mask);

mask2 = ((omega_array > 3.17) & (k_array < 0.38)) | ((k_array > 0.21) & (omega_array < 0.95));
mask2 = ~mask2;
omega_array = omega_array(mask2);
k_array = k_array(mask2);
A_array = A_array(mask2);
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
hold on 
loglog(k_list,harmonic3,'--g')

axis([0.05 4 , 1e-1 30])
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
closest_harmonic(k_array < 0.41) = 1;

% set lowest values of harmonic2 to closest_harmonic = 2
% figure,
% mask = (k_array < 0.62) & (closest_harmonic == 3);
% loglog(k_array(mask),omega_array(mask),'o')
% axis([1e-1 2, 5e-1 5])
closest_harmonic((k_array < 0.62) & (closest_harmonic == 3)) = 2;

%% Plot harmonics with different colours 

% colors = ["blue","red","green"];
colors = ["#0072BD","#D95319","#77AC30"];

figure(15)
for i = 1 :3
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
hold on 
loglog(k_list,harmonic3,'--g')
axis([0.05 4 , 2e-1 10])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;
% 
% figname = [fig_folder 'Harmonics_shallow_water_hw_' replace(num2str(h_w),'.','p') '_different_colors'];
% saveas(gcf,figname,'fig')
% saveas(gcf,figname,'pdf')
% saveas(gcf,figname,'png')

%% Save selected points

filename = 'Data_plot_selected_harmonics_0226_mesange_waves_005';
directory = 'W:/SagWin2024/Data/0226/Drones/mesange/matData/10-waves_005/Plots/';
save([directory filename],'omega_array','k_array','A_array','closest_harmonic','colors','k_list','harmonic1','harmonic2','harmonic3','-v7.3')

%% Move each branch back to main one 

figure(16)
for i = 1 :3
    
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
hold on 
loglog(k_list,harmonic3,'--g')
axis([0.05 4 , 2e-1 10])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;

% figname = [fig_folder 'Recomposition_harmonics_water_hw_' replace(num2str(h_w),'.','p')];
% saveas(gcf,figname,'fig')
% saveas(gcf,figname,'pdf')
% saveas(gcf,figname,'png')


% #######################################
%% Attenuation coefficient in corrected direction 
% #######################################

% Parameters 
selected_freq = [0.15 0.33]; % selected frequencies between which we proceed to the analysis
% get indices of frequencies closest to max and min selected frequency
[min_freq, start_freq_idx] = min(abs(f - selected_freq(1))); 
[max_freq, end_freq_idx] = min(abs(f - selected_freq(2)));
new_freq = f(start_freq_idx : end_freq_idx); % new frequency array 
FFT_cropped = FFT_t(:,:,start_freq_idx:end_freq_idx);
nf = length(new_freq);

x = m.x;
y = m.y;
padding_bool = 1;
add_pow2 = 2;
threshold = 0.8;
L0 = 130; % segments size in meter 

new_folder_fig = [fig_folder 'attenuation_oriented/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

alpha = zeros(nf,1); % array of attenuation coefficients
d = zeros(nf,1); % array of distance to plot

for idx_freq = 1:nf

    field = FFT_cropped(:,:,idx_freq);
    freq_txt = replace(num2str(new_freq(idx_freq)),'.','p');
    % 2D FFT of this real field 
    [shifted_fft,fft_2D,kx,ky] = spatial_FFT(field,padding_bool,add_pow2,facq_x);

    % Detect peaks 
    zmax = max(abs(shifted_fft),[],'all');
    norm_fft = abs(shifted_fft)/zmax;

    i = 1; % select first peak 
    % binarize the FFT spectrum
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

    % Build several segments for theta < 0
    % initial line 
    x0 = m.x(1); % minimal value along x axis
    y0 = L0*sin(abs(theta)) + m.y(end);
    ds = m.fx; % step of curvilinear coordinate

    s = (0:ds:L0);
    % Points that define the line 
    x_line = x0+s*cos(theta);
    y_line = y0+s*sin(theta);

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
    
    % Compute field correctly oriented
    [field_star,x_star,y_star] = orientation_parallel2propagation(field,theta,m.fx,L0);
    
    % Exponential fit 
    TF_x = squeeze(mean(field_star,2)); % average the TF over y  
    A = abs(TF_x); % amplitude along x-axis for each frequency
    
    % indices used to fit the exponential decay
    i_min = 1;
    i_max = length(A);

    decay_fig = figure(2);
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

%% Save relevant variables 
attenuation_file =  ['Data_attenuation_oriented_' num2str(selected_freq(1)) 'Hz_to_' num2str(selected_freq(2)) 'Hz'];
attenuation_file = replace(attenuation_file,'.','p');
attenuation_file = [fig_folder attenuation_file '_mode_1'];
save(attenuation_file,'alpha','d','new_freq','selected_freq','L0','threshold','padding_bool','add_pow2')

disp('DONE.')

% ##################################
%% Attenuation coefficient mode 1 
% ##################################

% load parameters 
base = 'W:/SagWin2024/Data/0226/Drones/mesange/matData/10-waves_005/Plots/';
filename = [base 'attenuation_coef_data_fmin0p1_fmax0p8_fthresh0p36 _mode1.mat'];
S1 = load(filename);

filename = [base 'Data_attenuation_oriented_0p15Hz_to_0p33Hz_mode_1.mat'];
S2 = load(filename);

%% merge data sets 
fmin = 0.15; 
fmax = 0.33; 

freq_array = zeros(size(S1.freq));
alpha = zeros(size(S1.lambda));
dist = zeros(size(S1.dist_fit));

for i = 1:length(freq_array)
    current_f = S1.freq(i);
    if ((S1.freq(i) >= fmin) && (S1.freq(i) <= fmax))
        [~,i0] = min(abs(S2.new_freq - S1.freq(i)));
        freq_array(i) = S2.new_freq(i0);
        alpha(i) = S2.alpha(i0);
        dist(i) = S2.d(i0);
    else
        freq_array(i) = S1.freq(i);
        alpha(i) = S1.lambda(i);
        dist(i) = S1.dist_fit(i);
    end 
    
end 

freq_thresh = S1.freq_thresh;
x_bound_thresh = S1.x_bound_thresh;
savename = [base 'Data_attenuation_merged_0p1Hz_to_0p8Hz_mode_1.mat'];
save(savename,'freq_array','alpha','dist','freq_thresh','x_bound_thresh')
%%

%  masking according to dist_fit
d_thresh = 0.05;

mask = (dist < d_thresh) & (abs(alpha) > 0.001) ; % keep only points for which distance is smaller than..
f = freq_array;
fitted_f = f(mask);
fitted_alpha = abs(alpha(mask))';

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
attenuation_filename = [fig_folder 'attenuation_law_mode_1_merged_confidence_' thresh_txt];
saveas(attenuation_fig,attenuation_filename,'fig');
saveas(attenuation_fig,attenuation_filename,'pdf');


% ####################
%% MAIN DEVELOPMENTS 
% ####################

% ######################################################
%% Try to fit a lorentzian on each peak of the 2D - FFT 
% ######################################################

% select a frequency 
selected_freq = 0.5;
[~,idx_freq] = min(abs(f - selected_freq));

field = FFT_t(:,:,idx_freq);

figure(17)
x = m.x;
y = m.y;

pcolor(x,y,real(field)')
shading interp
colormap(redblue)
axis image

% plot a profile 
selected_y = 50; % in meter
[~,i_y] = min(abs(y - selected_y));
figure(18)
plot(x,real(field(:,i_y)),'.')

%% 2D FFT of this real field 

padding_bool = 1;
add_pow2 = 2;
[shifted_fft,fft_2D,kx,ky] = spatial_FFT(field,padding_bool,add_pow2,facq_x);

figure(19)
pcolor(kx,ky,abs(shifted_fft)')
shading interp

[KY,KX] = meshgrid(ky,kx);
figure(20)
surf(KX,KY,abs(shifted_fft)')
shading interp
axis([-2 2, -2 2])
set(gca,'ColorScale','log')
%% Detect peaks 

zmax = max(abs(shifted_fft),[],'all');
norm_fft = abs(shifted_fft)/zmax;

figure(21)
surf(KX,KY,norm_fft')
shading interp
axis([-2 2, -2 2])
set(gca,'ColorScale','log')
%%
% binarize the FFT spectrum
threshold = 0.6;
bin_fft = norm_fft;
bin_fft(norm_fft > threshold) = 1;
bin_fft(norm_fft <= threshold) = 0;

CC = bwconncomp(bin_fft');
% stats = regionprops("table",CC,'Centroid');

stats = regionprops("table",CC,norm_fft','WeightedCentroid','MaxIntensity','PixelIdxList');
% transposed = norm_fft';
% [~,idx_peak] = min(abs(transposed(stats.PixelIdxList{1}) - stats.MaxIntensity));
% [row,col] = ind2sub(size(transposed),stats.PixelIdxList{1}(idx_peak));
% kx_peak = kx(row);
% ky_peak = ky(col);

figure(22) 
imagesc(kx,ky,norm_fft')
hold on 
plot(kx(round(stats.WeightedCentroid(:,1))),ky(round(stats.WeightedCentroid(:,2))),'.r')
% plot(kx_peak,ky_peak,'.r')

%% Fit a Lorentzian curve on a peak 

[~,i0] = min(abs(freq - selected_freq));
disp(['f = ' num2str(freq(i0))])
disp(['From sinus fitting, alpha = ' num2str(lambda(i0))])

i = 1; % index of the detected peak 

row = round(stats.WeightedCentroid(i,1));
col = round(stats.WeightedCentroid(i,2));

% coordinates of the peak in Fourier space 
kx_peak = kx(row);
ky_peak = ky(col); 

% convert to polar coordinates 
[theta,k_peak] = cart2pol(kx_peak,ky_peak);

% draw a line passing through peak and origin 
L_seg = 0.15; % half size of the segment, in k units
kx0 = kx_peak - cos(theta)*L_seg;
ky0 = ky_peak - sin(theta)*L_seg; 
ds = 2*pi*facq_x/size(shifted_fft,1); % step of the curvilinear coordinate 
s = (0:ds:2*L_seg); % curvilinear coordinate 
kx_line = kx0 + s*cos(theta);
ky_line = ky0 + s*sin(theta); 

figure(23),
pcolor(kx,ky,abs(shifted_fft)')
shading interp
hold on 
plot(kx_line,ky_line,'r--')
hold on 
plot(kx_peak,ky_peak,'r.')
hold on 
plot(kx(size(shifted_fft,1)/2 + 1),ky(size(shifted_fft,2)/2 + 1),'g.')
axis image 

%% Interpolate Fourier spectrum on this line
kx_inv = sort(kx); % sort kx in increasing order 
ky_inv = sort(ky);
F = griddedInterpolant({kx_inv,ky_inv},norm_fft); % create interpolant 

values = F(kx(end) - kx_line + kx(1),ky_line);
if ~rem(length(s),2)
    middle_s = length(s)/2 + 1;
else
    middle_s = (length(s)+1)/2 ;
end 
s = s - s(middle_s);
% fit by a lorentzian 
[y_lorentz_0,p0,~] = lorentzfit(s,values,[],[],'2');
figure(24)
plot(s,abs(values),'.')
hold on 
plot(s(middle_s),abs(values(middle_s)),'dr')
hold on 
plot(s,y_lorentz_0,'r--')

alpha = sqrt(p0(2));
disp(['Alpha = ' num2str(alpha) ' m-1'])
%%
% select a profile along kx 

window = 15;
profile_kx = abs(shifted_fft(row - window : row + window,col)).^2;
kx_tofit = kx(row - window : row + window);
kx_tofit = kx_tofit - kx(row);

% kx_inv = kx_tofit(end:-1:1);
% profile = profile_kx(end:-1:1);
% % Interpolate data
% F = griddedInterpolant(kx_inv,profile);
% k_precise = linspace(0.15,0.45,100);
% value = F(k_precise);

% use lorentzfit function 
P0 = [abs(shifted_fft(row))^2, kx(row), 0.1]; % starting values of fitting parameters
[y_lorentz_0,p0,~] = lorentzfit(kx_tofit,profile_kx',[],[],'2');

% use home made lorentzfit function
P1 = fminsearch(@(s)lorentzian_fit_2param(kx_tofit,profile_kx',s),[1,1]);

k_list = linspace(kx_tofit(1),kx_tofit(end),100);
yth = lorentzian_fun_2param(k_list,P1); % fitted exponential function

figure(24),
plot(kx_tofit,profile_kx,'.')
hold on 
plot(kx_tofit,y_lorentz_0,'r-')
hold on 
plot(k_list,yth,'b--')

alpha = sqrt(p0(2));
disp(['Alpha = ' num2str(alpha)])

alpha_2 = sqrt(P1(2));
disp(['Alpha = ' num2str(alpha_2)])


%%
figure,
pcolor(kx,ky,abs(A)')
shading interp
xlabel('$k_x \: \rm (m^{-1})$')
ylabel('$k_y \: \rm (m^{-1})$')
colorbar()
hold on 
plot(kx(x0),ky(y0),'dr')

figure,
pcolor(abs(A)')
% shading interp
xlabel('$k_x \: \rm (m^{-1})$')
ylabel('$k_y \: \rm (m^{-1})$')
colorbar()
hold on 
plot(x0,y0,'dr')

%% Plot radii over which we performed velocity field 

R = 2*pi*R_tics*facq_x/padding(1);
theta = linspace(-pi,pi,50);
kx_rad = R.*cos(theta);
ky_rad = R.*sin(theta);

%% 
figure, 
plot(R_tics,R_profile)
grid on 

%% 
[~,idx] = max(R_profile);

figure, 
pcolor(kx,ky,abs(A)')
% shading interp
xlabel('$k_x \: \rm (m^{-1})$')
ylabel('$k_y \: \rm (m^{-1})$')
colorbar()
axis image
for idx = 1:2
    hold on 
    plot(kx_rad(idx,:),ky_rad(idx,:),'r.')
end 


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
