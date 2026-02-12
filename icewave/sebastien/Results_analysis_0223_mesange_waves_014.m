clear all
close all

%% Add path to matlab functions 

functions_path = 'C:\Users\sebas\icewave\icewave\drone\Drone_banquise_analysis\';
% Add the folder and all its subfolders to the MATLAB path
addpath(genpath(functions_path));

%% Load structured data 
year = '2024';
date = '0223';
drone = 'mesange';
exp_ID = '35-waves_014';
suffixe_fig = [year '_' date '_' drone '_' exp_ID] ;


path2data = ['H:/Rimouski_2024/Data/' year '/' date '/' drone '/matData/' exp_ID '/'];
filename = 'PIV_processed_i00_Dt4_b1_W32_xROI650_width3190_yROI1_height_2159_scaled.mat';
disp('Loading Data..')
load([path2data filename])
disp('Data loaded')

%% Create a folder where figures are saved and define parameters
fig_folder = [path2data 'Plots/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end

font_size = 13;
facq_x = 1/m.SCALE.fx; % scale in box / meter
facq_t = 1/m.SCALE.ft; % scale in frame / sec
scale_V = m.SCALE.scale_V; % factor scaling for velocity in meter / s
W = m.PIV_param.w ;

%% Get Histogram displacement 

filename = [fig_folder 'Histogramm_displacements_Vx_time_average_' year '_' date '_' drone '_' exp_ID];
average_bool = 1;
get_histogram_displacement(m.Vx/scale_V,W,average_bool,font_size,filename);

%% Get profile of velocity
i_x = 30;
i_y = 60;
left_bool = 1;
% Vx = flip(m.Vx,1);
profile_fig = plot_located_profile(m.Vx,i_x,i_y,facq_x,facq_t,left_bool);
set_Papermode(profile_fig)

fig_file_name = [fig_folder 'Time_profile_'   '_ix_' num2str(i_x) '_iy_' num2str(i_y)];
saveas(profile_fig,fig_file_name,'fig')

%% Get rid off quadratic noise (drone movements)
[nx,ny,nt] = size(m.Vx);
x = (1:1:nx);
y = (ny:-1:1);

% compute the mean component of the velocity field for each frame
Vxmoy = mean(mean(m.Vx,2),1);
Vymoy = mean(mean(m.Vy,2),1);

Vx = supress_quadratic_noise(m.Vx,x,y); 
Vy = supress_quadratic_noise(m.Vy,x,y);
% Vx = flip(Vx,1); % waves are initially coming from the right boarder of the image
% Vy = flip(Vy,1);

disp('Drone motion corrected')

%% Plot main features of the velocity field
disp('Get velocity features')
idx_frame = 20;
caxis_amp = 1.0; % Amplitude of the colorbar scale (in meter / second)

plot_velocity_features(Vx,Vy,facq_x,idx_frame,caxis_amp,fig_folder)

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

file = [fig_folder 'FFT_spectrum_' suffixe_fig];
saveas(fig_spectrum,file,'fig')
saveas(fig_spectrum,file,'pdf')

%% Get demodulated field
disp('Getting demodulated fields')
selected_freq = [0.1 1.0];
x_bound = [1 size(FFT_t,1)];
caxis_amp = 0.2; % amplitude of the colorbar in meter/second
fig_name = ['Demodulated_field_Vx_' suffixe_fig 'caxis_' num2str(caxis_amp) 'ms'];

save_image = 1;
save_video = 1;
plot_demodulated_field(FFT_t,f,facq_x,selected_freq,x_bound,caxis_amp,fig_folder,fig_name,save_image,save_video)


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
disp('E(k,f) computed')

k = 2*pi*R_tics*facq_x/padding(1);

%% Plot A(k,f)
omega = 2*pi*f;
k_list = linspace(0.1,7,100);

g = 9.81; % gravity intensity 
h_w = 4.8; % water depth
rho_ice = 917; % ice density in kg.m-3
rho_w = 1000; % water density in kg.m-3
h_ice = 0.2; % ice thickness in meter

deep_water = sqrt(g*k_list);
shallow = sqrt(g*h_w*k_list.^2);
yth = sqrt(g*k_list.*tanh(h_w*k_list));

y_squire = sqrt(g*k_list.*tanh(h_w*k_list)./(1 + (rho_ice/rho_w)*k_list*h_ice.*tanh(k_list*h_w)));

harmonic2 = sqrt(2*g*k_list.*tanh(h_w*k_list/2));
harmonic3 = sqrt(3*g*k_list.*tanh(h_w*k_list/3));

% plot A(omega,k)
fig_Afk = figure; 
pcolor(k,omega,E)
shading interp
xlabel('$k \: \rm (rad.m^{-1})$')
ylabel('$\omega \: \rm (rad.s^{-1})$')
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'ColorScale','log')
% colormap(slanCM('gnuplot2'))
cbar = colorbar();
cbar.Label.String = '$|\hat{V}_x|(k,\omega)$';
cbar.Label.Interpreter = 'latex';
axis([0.1 5, 0.7 7.0])
hold on 
loglog(k_list,yth,'w--')
hold on
loglog(k_list,harmonic2,'w--')
% hold on 
% loglog(k_list,y_squire,'r--')
hold on 
loglog(k_list,harmonic3,'w--')
set_Papermode(gcf)
set(gca,'FontSize',13)
clim([1e-4 2e-2])
% xticks([1e-1 2e-1 3e-1 1 2 3])
% lgnd = legend('',['$\omega^2 = gk \tanh(gh_w) \: h_w = ' num2str(h_w) '\: \rm m$'],'','');
% set(lgnd,'Location','southeast')
% set(lgnd,'color','none')
% loglog(k_list,harmonic2,'w--')
% hold on 
% loglog(k_list,harmonic3,'w--')
% colormap(slanCM('thermal-2'))

fig_filename = ['Afk_' suffixe_fig 'with_harmonics_hw' num2str(h_w)];
fig_filename = replace(fig_filename,'.','p');

saveas(fig_Afk,[fig_folder fig_filename],'fig')
saveas(fig_Afk,[fig_folder fig_filename],'pdf')

%% Save data
data_filename = ['Data_Afk_' suffixe_fig];
data_filename = replace(data_filename,'.','p');

save([fig_folder data_filename],'f','omega','k','E','shift','-v7.3')
disp('Data saved')


%% Load Data for A(f,k) plot
data_filename = ['Data_Afk_' suffixe_fig];
data_filename = replace(data_filename,'.','p');

disp('Loading data..')
load([fig_folder data_filename])
disp('Data loaded')


%% Detect each branch of bound waves

minimal_prominence = 1.5e-1;
freq_end = 1.0; % stop loop at this frequency
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
    [pks,locs,w,prom] = findpeaks(profile,'MinPeakProminence',minimal_prominence);
    [y_max,k_peaks] = subpix_precision(profile,k',locs);
    
    M_peaks(i).k = k_peaks;
    M_peaks(i).A = y_max*max_absolute;
    M_peaks(i).width = w;
    M_peaks(i).omega = omega(i);

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

mask = ~(((omega_array > 2.2) & (k_array < 0.45)) | (omega_array < 0.55)); 
omega_array = omega_array(mask);
k_array = k_array(mask);
A_array = amplitude_array(mask);


mask2 = ((omega_array > 1.02) & (k_array < 0.15)) | ((k_array > 0.35) & (omega_array < 1.3))...
    | ((omega_array < 0.575) | ((k_array < 0.171) & omega_array > 1.195)) | ((omega_array > 3.98) & (k_array < 0.68)) ...
    | k_array < 0.1 | (k_array < 0.55) & (omega_array > 3.0);
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
h_w = 4.8; % water depth
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

% axis([0.05 4 , 1e-1 30])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;


%% Associate each points to an harmonic 

% compute a distance to each harmonics 
N = [1,2,3]';
% k = [1 1.5 2];

% #idx1 : N, idx2 : k
omegaN = bound_harmonicN(k_array,N,h_w);
D = sqrt((repmat(omega_array,[3,1]) - omegaN).*2); 
[~,closest_harmonic] = min(D,[],1);

% set lowest values of k to closest_harmonic = 1
closest_harmonic(k_array < 0.45) = 1;

% set lowest values of harmonic2 to closest_harmonic = 2
figure,
mask = (k_array < 0.6) & (closest_harmonic == 3);
loglog(k_array(mask),omega_array(mask),'o')
axis([1e-1 2, 5e-1 5])
closest_harmonic((k_array < 0.62) & (closest_harmonic == 3)) = 2;

%% Plot harmonics with different colours 

colors = ["#0072BD","#D95319","#77AC30"];
% colors = ["#0072BD","#D95319"];
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
% axis([0.05 4 , 2e-1 10])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;


%% Save selected points

file2save = [ fig_folder 'Data_plot_selected_harmonics_' suffixe_fig];
save(file2save,'omega_array','k_array','A_array','closest_harmonic','colors','k_list','harmonic1','harmonic2', ...
    'harmonic3','-v7.3')

%% Plot each harmonic with a different colormaps corresponding to the intensity of each points (omega,k)

% Load Data of selected harmonics 
file2load = [ fig_folder 'Data_plot_selected_harmonics_' suffixe_fig];
load(file2load)
disp('Data of separated harmonics loaded')

%%
A_h1 = A_array(closest_harmonic == 1);
A_h1 = (A_h1 - min(A_h1))/(max(A_h1) - min(A_h1));

A_h2 = A_array(closest_harmonic == 2);
A_h2 = (A_h2 - min(A_h2))/(max(A_h2) - min(A_h2));

A_h3 = A_array(closest_harmonic == 3);
A_h3 = (A_h3 - min(A_h3))/(max(A_h3) - min(A_h3));
marker_size = 50;

harmonic_cmaps = figure;
% create two different axes 
ax1 = axes(harmonic_cmaps);
ax2 = copyobj(ax1,harmonic_cmaps);
ax3 = copyobj(ax1,harmonic_cmaps);

plot(ax1,k_list,harmonic1,'--b')
hold on
plot(ax1,k_list,harmonic2,'--r')
hold on 
plot(ax1,k_list, harmonic3,'--g')
s1 = scatter(ax1,k_array(closest_harmonic == 1),omega_array(closest_harmonic == 1),marker_size,A_h1,'filled','MarkerEdgeColor','k');

s2 = scatter(ax2,k_array(closest_harmonic == 2),omega_array(closest_harmonic == 2),marker_size,A_h2,'filled','MarkerEdgeColor','k');

s3 = scatter(ax3,k_array(closest_harmonic ==3), omega_array(closest_harmonic == 3), marker_size, A_h3, 'filled', ...
    'MarkerEdgeColor','k');
colormap(ax1,slanCM(9))
colormap(ax2,slanCM(11))
colormap(ax3,slanCM(10))
ax2.Visible = "off";
ax3.Visible = "off";

% sets limits of each axes
xmin = 6e-2;
xmax = 5;
ymin = 0.3;
ymax = 8;
xlim(ax1,[xmin xmax])
xlim(ax2,[xmin xmax])
xlim(ax3,[xmin xmax])
ylim(ax1,[ymin ymax])
ylim(ax2,[ymin ymax])
ylim(ax3,[ymin ymax])
grid on 
% set axes scale in log-log
set(ax1,'xscale','log')
set(ax1,'yscale','log')
set(ax2,'xscale','log')
set(ax2,'yscale','log')
set(ax3,'xscale','log')
set(ax3,'yscale','log')

xlabel(ax1,'$k \: \rm (rad.m^{-1})$')
ylabel(ax1,'$\omega \: \rm (rad.s^{-1})$')
legend(ax1,['$\omega_1 = \sqrt{gk \tanh(' num2str(h_w) 'k)}$'],...
    ['$\omega_2 = \sqrt{2gk \tanh(' num2str(h_w) '\frac{k}{2})}$'],...,
    ['$\omega_3 = \sqrt{3gk \tanh(' num2str(h_w) '\frac{k}{3})}$'],'','Location','southeast')
ax1.FontSize = 13;
set_Papermode(harmonic_cmaps)

fig_filename = [fig_folder 'Harmonics_colormaps_' suffixe_fig '_hw_' replace(num2str(h_w),'.','p') ];
saveas(harmonic_cmaps,fig_filename,'fig')
saveas(harmonic_cmaps,fig_filename,'pdf')

%% Move each branch back to main one 
recompo_harmonic = figure;
% create two different axes 
ax1 = axes(recompo_harmonic);
ax2 = copyobj(ax1,recompo_harmonic);
ax3 = copyobj(ax1,recompo_harmonic);

plot(ax1,k_list,harmonic1,'--b')
hold on
plot(ax1,k_list,harmonic2,'--r')
hold on 
plot(ax1,k_list, harmonic3,'--g')
s1 = scatter(ax1,k_array(closest_harmonic == 1),omega_array(closest_harmonic == 1),marker_size,A_h1,'filled','MarkerEdgeColor','k');

s2 = scatter(ax2,k_array(closest_harmonic == 2)/2,omega_array(closest_harmonic == 2)/2,marker_size,A_h2,'filled','MarkerEdgeColor','k');

s3 = scatter(ax3,k_array(closest_harmonic ==3)/3, omega_array(closest_harmonic == 3)/3, marker_size, A_h3, 'filled', ...
    'MarkerEdgeColor','k');
colormap(ax1,slanCM(9))
colormap(ax2,slanCM(11))
colormap(ax3,slanCM(10))
ax2.Visible = "off";
ax3.Visible = "off";

% sets limits of each axes
xmin = 6e-2;
xmax = 5;
ymin = 0.3;
ymax = 8;
xlim(ax1,[xmin xmax])
xlim(ax2,[xmin xmax])
xlim(ax3,[xmin xmax])
ylim(ax1,[ymin ymax])
ylim(ax2,[ymin ymax])
ylim(ax3,[ymin ymax])
grid on 
% set axes scale in log-log
set(ax1,'xscale','log')
set(ax1,'yscale','log')
set(ax2,'xscale','log')
set(ax2,'yscale','log')
set(ax3,'xscale','log')
set(ax3,'yscale','log')


set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;
legend(ax1,['$\omega_1 = \sqrt{gk \tanh(' num2str(h_w) 'k)}$'],...
    ['$\omega_2 = \sqrt{2gk \tanh(' num2str(h_w) '\frac{k}{2})}$'],...,
    ['$\omega_3 = \sqrt{3gk \tanh(' num2str(h_w) '\frac{k}{3})}$'],'','Location','southeast')
figname = [fig_folder 'Recomposition_harmonics_water_' suffixe_fig '_hw_' replace(num2str(h_w),'.','p')];
saveas(gcf,figname,'fig')
saveas(gcf,figname,'pdf')
saveas(gcf,figname,'png')



%%

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
axis([xmin xmax ymin ymax])


% #######################################
%% Attenuation coefficient in corrected direction 
% #######################################

% Parameters 
selected_freq = [0.15 0.5]; % selected frequencies between which we proceed to the analysis
% get indices of frequencies closest to max and min selected frequency
[min_freq, start_freq_idx] = min(abs(f - selected_freq(1))); 
[max_freq, end_freq_idx] = min(abs(f - selected_freq(2)));
new_freq = f(start_freq_idx : end_freq_idx); % new frequency array 
FFT_cropped = FFT_t(:,:,start_freq_idx:end_freq_idx);
nf = length(new_freq);

% Select range of frequencies and space of shortened fit 
f_short = [0.45]; % frequency above which the attenuation fit is shortened
x_short = [80]; % meters % range over which the shortened attenuation fit is performed

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
C  =zeros(nx,1);
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
    C(idx_freq) = 10.^p(0); % prefactor

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
attenuation_file =  ['Data_attenuation_oriented' suffixe_fig '_' num2str(selected_freq(1)) 'Hz_to_' ...
    num2str(selected_freq(2)) 'Hz' '_mode_1_analyse_1021'];
attenuation_file = replace(attenuation_file,'.','p');
attenuation_file = [fig_folder attenuation_file];
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
d_thresh = 0.3;

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
attenuation_filename = [fig_folder 'attenuation_law_ ' suffixe_fig '_mode_1_confidence_' thresh_txt];
saveas(attenuation_fig,attenuation_filename,'fig');
saveas(attenuation_fig,attenuation_filename,'pdf');


%% Plot distance d to the curve as function of frequencies
figure, 

plot(S.new_freq,S.d,'o')
xlabel('$f \: \rm (Hz)$')
ylabel('$d$')




% ##########################
%% ######## FUNCTION SECTION 
% ##########################

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