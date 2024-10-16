

clear all; 
close all;

%% Loading structure obtained after PIV processing and post-processing

% base = 'E:/Rimouski_2024/Data/2024/0211/Drones/Fulmar/FULMAR_vertical/matdata/';
base = 'W:/SagWin2024/Data/0211/Drones/Fulmar/matData/';
%base = 'E:/PIVlab_drone/matdata/DJI_0308_Dt4_W64_full/';

filename = 'PIV_processed_i04500_Dt4_b1_W32_full_total_processed_scaled.mat';
matname = [base filename];
%%
disp('Loading Data..');
load(matname);
disp('Data loaded');
%% Creates a folder where to save the generated plots
base_fig = base;
%base_fig = 'E:/PIVlab_drone/';
fig_folder = [base_fig 'Plots/continuous_ice/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%%
w = m.w ;
L_x = 3840; % size of the image in pixel, along larger axis (x-axis)
h_drone = 140; % height of the drone in meter
theta_x = 32.75; % semi AFOV of the drone, along x-axis, in °

facq_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
facq_x = facq_pix*2/w; % scale in box / meter
fx = 1/facq_x; % factor scaling in meter / box
facq_t = 29.97; 

%% get histogramm

filename = [fig_folder 'histogramm_displacements_Vy_ice'];
font_size = 13;
bool_average = 0;

get_histogram_displacement(m.Vy(1:38*facq_x,:,:)/m.scale_V,w,bool_average,font_size,filename);

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
V = sqrt(Vx.^2 + Vy.^2);
disp('Drone motion corrected')

%% Get profile of velocity
i_x = 20;
i_y = 20;
plot_located_profile(V,i_x,i_y,facq_x,facq_t,fig_folder);

disp(mean(abs(V(i_x,i_y,:))))

%% Get time Fourier transform
disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;

[FFT_t,TF_spectrum,f] = temporal_FFT(V(:,:,:),padding_bool,add_pow2,facq_t);

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

%% Select a demodulated field for 2D-FFT

selected_freq = 0.39;
[val,idx] = min(abs(f - selected_freq));
% idx = idx + 5;
disp(f(idx))

xmin = 1;
xmax = 38*facq_x ;
a = FFT_t(xmin:xmax,:,idx);

padding_bool = 1;
add_pow2 = 2;

figure, 
imagesc(real(a)')
axis image

figure, 
imagesc(real(FFT_t(:,:,idx))')
axis image

% 2D-FFT 
[shifted_fft,fft_2D,kx,ky] = spatial_FFT(a,padding_bool,add_pow2,facq_x);
shifted_fft = shifted_fft ./ max(shifted_fft,[],'all');

mask_x = 2;
mask_y = 12;
shifted_fft(size(shifted_fft,1)/2+1-mask_x:size(shifted_fft,1)/2+1+mask_x,size(shifted_fft,2)/2+1-mask_y:size(shifted_fft,2)/2+1+mask_y) = zeros;

figure, imagesc(kx,ky,abs(shifted_fft)')

% Detect peaks 
min_height = 0.8;
[pks,locs_x,locs_y] = peaks2(shifted_fft,'MinPeakHeight',min_height);
kx_peaks = kx(locs_x);
ky_peaks = ky(locs_y);
%%
% Plot detected peaks
figure,
imagesc(kx,ky,abs(shifted_fft)');
shading interp
hold on 
plot(kx_peaks,ky_peaks,'ro')















%% Get demodulated field
disp('Getting demodulated fields')
selected_freq = [0.1 1.2];
x_bound = [1 38*facq_x];
caxis_amp = -1; % amplitude of the colorbar in meter/second
fig_name = 'Demodulated_field_continuous';

save_image = 0;
save_video = 0;
plot_demodulated_field(FFT_t,f,facq_x,selected_freq,x_bound,caxis_amp,fig_folder,fig_name,save_image,save_video)

%% Get wave vectors 
disp('Getting wave vectors')
selected_freq = [0.1 1.2]; % selected frequencies between which we proceed to the analysis
x_bound = [1 38*facq_x]; % selected boundaries at which we perform 2D FFT
padding_bool = 1;
add_pow2 = 2; % additional power of 2 for padding 
black_mask = 7;
caxis_amp = -1;
fig_name = 'Spatial_Fourier_space_video_V';

save_image = 0;
save_video = 0;

[k,freq] = get_wave_vectors(FFT_t,f,facq_x,selected_freq,x_bound,padding_bool,add_pow2,black_mask,caxis_amp,fig_folder,fig_name,save_image,save_video);

% save data to plot dispersion relation
filename = ['dispersion_relation_data_continuous_fmin' num2str(selected_freq(1)) '_fmax' num2str(selected_freq(2)) '_add_pow' num2str(add_pow2)];
filename = replace(filename,'.','p');
dispersion_file = [fig_folder filename];
save(dispersion_file,'k','freq','selected_freq','x_bound','add_pow2','black_mask')

%%

omega = 2*pi*freq;
wave_vector = k;

% first plot before masking
figure,
loglog(wave_vector,omega,'o');
grid on 
xlabel('$k$ $(m^{-1})$','Interpreter','latex');
ylabel('$\omega$ $(s^{-1})$','Interpreter','latex');

%% Computing filtered datas
% masking some datas
%mask = (wave_vector > 0.12) & (omega < 2.78);
opposite_mask = (omega > 3.5) & (wave_vector < 0.2);
mask = ~opposite_mask;
fitted_omega = omega(mask);
fitted_k = wave_vector(mask);

figure,
loglog(fitted_k,fitted_omega,'o');
grid on 
xlabel('$k$ $(m^{-1})$','Interpreter','latex');
ylabel('$\omega$ $(s^{-1})$','Interpreter','latex');

%% Final Plot !

data_plot = [fig_folder 'Disp_relation_addpad2.mat'];
load_bool = 0;
if load_bool
    load(data_plot)
    disp('Data to plot loaded');
end 

save_boolean = 0;

%Physical_parameters
g = 9.81; % intensity of gravity
h = 1.00; % water depth in meter
E = 2.5e9;
h_ice = 12*0.01; 
nu = 0.28;
rho_w = 1000;
D = E*h_ice^3/12/rho_w/(1-nu^2); % Flexion modulus
rho = 1000; % water density

% minimal values of the plot on x and y axis
xlim_min = 0.1;
xlim_max = 10;
ylim_min = 0.1;
ylim_max = 20;

% Computation of the different models 
k_list = linspace(xlim_min,xlim_max,1000); % array of k for plots
deep_water = sqrt(g*k_list);
shallow_water = sqrt(g*h*k_list.^2);
flexural = sqrt(g.*k_list + D.*(k_list.^5));
gravito_waves = sqrt(g*k_list .* tanh(k_list*h));

% Power law fitting 
% data to plot
mask_powerlaw = (fitted_omega < 2.0) & (fitted_k > 0.53) ;
k_powerlaw = fitted_k(mask_powerlaw);
omega_powerlaw = fitted_omega(mask_powerlaw);
l1 = fminsearch(@(s)powerfit(k_powerlaw,omega_powerlaw,s),[1,1]);
% l1 = fminsearch(@(s)powerfit(fitted_k,fitted_omega,s),[1,1]);
yth = powerfun(k_list,l1); % fitted exponential function

% Plot
relation_disp_scaled = figure(34);
%loglog(wave_vector,new_omega,'o');
loglog(fitted_k,fitted_omega,'o');
hold on 
loglog(k_list,deep_water,'k-');
hold on 
loglog(k_list,flexural,'b-');
% hold on 
% loglog(k_list,yth,'r-');
hold on 
loglog(k_list,gravito_waves,'k--');
% hold on 
% loglog(k_powerlaw,omega_powerlaw,'o','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor','black');
% title('Dispersion Relation');
xlabel('Wave number $k$ $(m^{-1})$','Interpreter','latex');
ylabel('Pulsation $\omega$ $(s^{-1})$','Interpreter','latex');
axis([xlim_min xlim_max ylim_min ylim_max]);
shallow_txt = ['Shallow, $h = ' sprintf('%0.1f',h) '\: \rm m$'];
power_law_txt = ['$\omega = ' sprintf('%0.2f',l1(2)) 'k^{' sprintf('%0.2f',l1(1)) '}$'];
gravito_txt = ['$\omega = \sqrt{(g  k  tanh(k h))}$'];
legend('Data','Deep water', shallow_txt, gravito_txt, 'Interpreter','latex','Location','northwest');
grid on

% Saving the plot
%fig_folder = 'C:/Users/sebas/Stage_MIZ/Traitement_PIV_Rimouski_20230310/';
file_relation_disp_scaled = [fig_folder 'Relation_disp_2024_02_11_ice'];
if save_boolean 
    saveas(relation_disp_scaled,file_relation_disp_scaled,'pdf');
    saveas(relation_disp_scaled,file_relation_disp_scaled,'fig');
end



%% FUNCTION SECTION

function yth = powerfun(x,l)
    yth = l(2)*x.^l(1);
end 

function d = powerfit(x,y,l)
    yth = powerfun(x,l);
    d = sum((yth-y).^2);
end 