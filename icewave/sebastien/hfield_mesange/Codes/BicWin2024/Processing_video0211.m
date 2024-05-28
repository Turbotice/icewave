%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

base = 'W:/SagWin2024/Data/0211/Drones/Fulmar/matData/';
%base = 'E:/PIVlab_drone/matdata/DJI_0308_Dt4_W64_full/';

filename = 'PIV_processed_i04500_Dt4_b1_W32_full_total_processed_Seb.mat';

% filename = 'PIV_processed_i04500_Dt4_b1_W32_full_total_processed.mat';
matname = [base filename];
%%
disp('Loading Data..');
load(matname);
disp('Data loaded');
%% Creates a folder where to save the generated plots
base_fig = 'W:/SagWin2024/Data/0211/Drones/Fulmar/matData/';
%base_fig = 'E:/PIVlab_drone/';
fig_folder = [base_fig 'SWO_FUL_20240211T203233UTC/Plots/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Scaling 
fe = 29.97; % Frame rate in Hz
nb_pass = m.s_param{6,2};
W = 32; % size of the window for DIC algorithm
font_size = 13;

% ##########################################
L_x = 3840; % size of the image in pixel, along x-axis
h_drone = 140; % height of the drone in meter
theta_x = 32.75; % semi AFOV of the drone, along x-axis, in °

fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
fx = fx_pix*2/W; % scale in box / meter 
% ##########################################

% fx = 0.8857; % for W = 64; 
%fx = 0.8857*2; % spatial scale in boxes/meter
Dt = 4; % step between two frames that were compared during the PIV algorithm 

scale_V = (fe/Dt) / fx_pix; % scale of the velocity in m/s

%% Save a scaled image 
picture_file = 'W:/SagWin2024/Data/0211/Drones/Fulmar/exemple_image.png';
img = imread(picture_file);

%%
% x = (1:1:size(img,2));
% y = (1:1:size(img,1));

font_size = 13;
new_img = img;
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

% figname = [fig_folder 'Im_exemple'];
% saveas(img_fig,figname,'fig')
% saveas(img_fig,figname,'pdf')
% saveas(img_fig,figname,'png')

%% Get Histogram displacement 

% Vx = m.Vx(10:end,:,:) - mean(mean(m.Vx(10:end,:,:),2),1);
filename = [fig_folder 'Histogram_displacement_Vy_continuous_ice_time_averaged'];
get_histogram_displacement(m.Vy(1:48,:,:)/m.scale_V,W,1,font_size,filename);

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

Vx = supress_quadratic_noise(m.Vx,x,y);
Vy = supress_quadratic_noise(m.Vy,x,y);


%% Show velocity field

for i=200:4:500
    figure(14)
    subplot(2,1,1)
    pcolor(x,y,m.Vy(:,:,i)')
    axis image 
    colormap(redblue)
    shading interp
    caxis([-1 1])
    colorbar()
    
    subplot(2,1,2)
    pcolor(x,y,m.Vy(:,:,i)')
    axis image 
    colormap(redblue)
    caxis([-1 1])
    shading interp
    colorbar()
    
    getframe();
    pause(0.2)
end

%% Get time Fourier transform
% We extract the demodulated wave field and try to project the velocity
% field on the orientation of propagation 
t = m.t;
x = m.x;
y = m.y;

% create a 3D grid [X,Y,T]
[X,Y,T] = ndgrid(x,y,t);

% sort the array of y scaled in ascending order
yi = sort(y);

% compute FFT in time 
disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;
[FFT_Vx,TF_spectrum_Vx,f] = temporal_FFT(Vx,padding_bool,add_pow2,fe);
[FFT_Vy,TF_spectrum_Vy,f] = temporal_FFT(Vy,padding_bool,add_pow2,fe);
disp('Time Fourier Transform computed')

%% Plot spectrum 

figure, 
loglog(f,TF_spectrum_Vx)
grid on 
xlabel('$f \: \rm (Hz)$')
ylabel('$\langle \hat {V_x} \rangle _{x,y}$');

figure, 
loglog(f,TF_spectrum_Vy)
grid on 
xlabel('$f \: \rm (Hz)$')
ylabel('$\langle \hat {V_y} \rangle _{x,y}$');

% ########################################
%% ####### Waves attenuation #############
% ########################################

% select demodulated wave field for a given frequency
selected_freq = 0.4; % selected frequencies between which we proceed to the analysis
[~, i0] = min(abs(f - selected_freq));
xmin = 1;
xmax = size(FFT_Vx,1);
disp(i0);
A_Vx = FFT_Vx(xmin:xmax,:,i0); % get the 2D map in temporal Fourier space
A_Vy = FFT_Vy(xmin:xmax,:,i0);

% computes wave vector components for a given frequency
padding_bool = 1;
add_pow2 = 2;
fx = 1/m.fx;
black_mask = 10;
caxis_amp = -1;
[~,kx_peak,ky_peak] = get_single_wavevector(A_Vx,padding_bool,add_pow2,fx,black_mask,caxis_amp);

[theta_Vx,k_Vx] = cart2pol(kx_peak,ky_peak); 

disp(['Theta = ' num2str(theta_Vx) ' rad'])
disp(['k = ' num2str(k_Vx) ' m-1'])

[~,kx_peak,ky_peak] = get_single_wavevector(A_Vy,padding_bool,add_pow2,fx,black_mask,caxis_amp);

[theta_Vy,k_Vy] = cart2pol(kx_peak,ky_peak); 

disp(['Theta = ' num2str(theta_Vy) ' rad'])
disp(['k = ' num2str(k_Vy) ' m-1'])

% final value of theta 
theta = mean([theta_Vx,theta_Vy]);

%% Show the demodulated field 
xmin = 1;
xmax = nx;%50*fx;

x = (xmin:1:xmax);
y = (ny:-1:1);
[X,Y] = meshgrid(x,y);

R = real(A_Vx(xmin:xmax,:)); % get real intensity map of FFT_temporal at a given frequency

figure,
pcolor(x/fx,y/fx,R')
% imagesc(x/fx,y/fx,R')
% 
% set(gca,'YDir','normal')
title(['Frequency : ' num2str(f(i0)) ' Hz'],'Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
shading interp
axis([0 size(X,2)/fx 0 size(X,1)/fx])
axis image
cbar = colorbar();
cbar.Label.String = '$ \rm{Re} \left( \overline {V_x}(x,y,f) \right) \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
%     cbar.Label.FontSize = font_size;
if caxis > 0 
    caxis([-caxis_amp caxis_amp])
end 
ax = gca;
ax.FontSize = 13;

%% Computes parallel and perpendicular velocities 
 
V_par = A_Vx*cos(theta) - A_Vy*sin(theta);
V_perp = A_Vx*sin(theta) + A_Vy*cos(theta);

caxis_amp = quantile(real(V_par),[0.05 0.95],'all');
% V_par = -A_Vx*cos(theta) - A_Vy*sin(theta);
% V_perp = A_Vx*sin(theta) - A_Vy*cos(theta);

figure, 
pcolor(m.x(xmin:xmax),m.y,real(V_par)')
shading interp
axis image 
colormap(redblue)
colorbar()
% caxis(caxis_amp)
xlabel('$x \: \rm (m)$','Interpreter','latex')
ylabel('$y \: \rm (m)$','Interpreter','latex')

figure,
pcolor(m.x(xmin:xmax),m.y,real(V_perp)')
shading interp
axis image 
colormap(redblue)
colorbar()
caxis(caxis_amp)
xlabel('$x \: \rm (m)$','Interpreter','latex')
ylabel('$y \: \rm (m)$','Interpreter','latex')

%% plot line in the direction of kx and ky

disp(['Theta = ' num2str(theta) ' rad'])

% coordinates of a point through which the line should pass (southeast
% point)

x0 = 60; % meter
y0 = 1; % meter
% y1 = abs(x0*tan(theta));
s = (0:+1:-x0/cos(theta)); % curvilinear coordinate

% Points that define the line 
x_line = x0+s*cos(theta);
y_line = y0+s*sin(theta);

figure(11)
hold off
pcolor(m.x,m.y,real(V_par)')
colormap(redblue)
% caxis(quantile(real(V_par),[0.05 0.95],'all'))
shading interp
hold on 
plot(x_line,y_line,'k--','LineWidth', 3)
axis image
% graphe_legende('$x$ (m)','$y$ (m)','$f$ = 0.3 Hz',true)

% saveas(gcf,[fig_folder 'champ_demodule_f20cHz_line'],'fig')
% saveas(gcf,[fig_folder 'champ_demodule_f20cHz_line'],'png')
% saveas(gcf,[fig_folder 'champ_demodule_f20cHz_line'],'eps')

%% Computes interpolant of parallel velocity
yi = sort(m.y);
Fpar = griddedInterpolant({m.x,yi,m.t},V_par);

% interpolate parallel velocity along the line
%% Interpolate pixels that are closest to this line

% y1 = y(y_line>y_scaled(end));
y1 = round(y_line);

%% Plot the projected velocity on the line 

x = x_scaled;
y = y_scaled;

% create a 3D grid [X,Y,T]
[X,Y,T] = ndgrid(x,y,t);

% sort the array of y scaled 
yi = sort(y);

% interpolation of Vx and Vy in the 3D space [x,y,t]
Fx = griddedInterpolant({x,yi,t},m.Vx);
Fy = griddedInterpolant({x,yi,t},m.Vy);


%test 
figure(12);
% coordinates (x0,y0) at which we display the projected velocity 
x0 = 60; % in meter 
y0 = 1/fx; % in meter 

% 
% [val,i0]= min(abs(x-x0)); % index for x-coordinate
% [val,j0]= min(abs(yi-y0)); % index for y-coordinate
% hold off
%plot(Fx(x0*ones(1,nt),y0*ones(1,nt),t))
%hold on
%plot(Fy(x0*ones(1,nt),y0*ones(1,nt),t))

Vx = Fx(x0*ones(1,nt),y0*ones(1,nt),t); % extract Vx using interpolated field 
Vy = Fy(x0*ones(1,nt),y0*ones(1,nt),t); % extract Vy using interpolated field 

% first test with Stephane
% plot(t,-Vx*cos(theta)-Vy*sin(theta))
% hold on
% plot(t,Vx*sin(theta)-Vy*cos(theta))
% hold on

%plot(t,Vx*cos(theta)-Vy*sin(theta))

% second test
% plot(t,+Vx*cos(theta) - Vy*sin(theta))
% hold on 
% plot(t,Vx*sin(theta) + Vy*cos(theta))


%%
for k=1:length(s)
    x0 = x_line(k);
    y0 = y_line(k);
    Lx(k,:)=Fx(x0*ones(1,nt),y0*ones(1,nt),t);
    Ly(k,:)=Fy(x0*ones(1,nt),y0*ones(1,nt),t);
    L_ll(k,:)=-Fx(x0*ones(1,nt),y0*ones(1,nt),t)*cos(theta);
end


%plot(Fx(100*ones(1,nt),40*ones(1,nt),t))

figure(13)
pcolor(t,-s,L1);shading interp
graphe_legende('$t$ s)','$L$ (m)','$f$ = 0.3 Hz',true);



% ####################################################
%% ########## Dispersion relation ####################
% ####################################################

% We extract the demodulated wave field and try to project the velocity
% field on the orientation of propagation 
t = (0:size(m.Vx,3))/fe;
x = m.x;
y = m.y;

% create a 3D grid [X,Y,T]
[X,Y,T] = ndgrid(x,y,t);

% sort the array of y scaled in ascending order
yi = sort(y);

% compute FFT in time 
disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;
[FFT_Vx,TF_spectrum,f] = temporal_FFT(Vx,padding_bool,add_pow2,fe);
[FFT_Vy,TF_spectrum,f] = temporal_FFT(Vy,padding_bool,add_pow2,fe);
disp('Time Fourier Transform computed')

%% Show the demodulated field 

fx = 1/m.fx; % scale in pix / meter
selected_freq = [0.1 0.8];
x_bound = [1 48]; % corresponds to ice 
caxis_amp = -1;
left_bool = 1;
fig_name = 'Demodulated_field_Vx_full';
save_image = 0;
save_video = 0;

plot_demodulated_field(FFT_Vy,f,fx,selected_freq,x_bound,caxis_amp,left_bool,fig_folder,fig_name,save_image,save_video)

%% Compute wave vectors for a given FFT 
padding_bool = 1;
add_pow2 = 2;
blac_mask = 10;
caxis_amp = -1;
save_image = 1;
save_video = 1;
vid_name = 'Continuous_ice_2D-FFT_Vx_addpow_2';
x_bound = [1 48]; % corresponds to ice (in boxes indices)

min_freq = 0.1; % minimal frequency in Hz
max_freq = 0.65; % maximal frequency in Hz 
selected_freq = [min_freq max_freq];
[dist,freq] = get_wave_vectors(FFT_Vx,f,fx,selected_freq,x_bound,padding_bool,add_pow2,black_mask,caxis_amp,fig_folder,vid_name,save_image,save_video);

%%
k = dist;
filename = ['dispersion_relation_data_continuous_fmin' num2str(selected_freq(1)) '_fmax' num2str(selected_freq(2)) '_add_pow' num2str(add_pow2) '_hamming'];
filename = replace(filename,'.','p');
dispersion_file = [fig_folder filename];
save(dispersion_file,'k','freq','selected_freq','x_bound','add_pow2','black_mask')


%% Plot dispersion relation 
% load data for dispersion relation 
base_disp = 'W:/SagWin2024/Data/0211/Drones/Fulmar/matData/SWO_FUL_20240211T203233UTC/Plots/';
filename_disp = 'dispersion_relation_data_continuous_fmin0p1_fmax0p65_add_pow2_hamming';

Sdisp = load([base_disp filename_disp]);
disp('Data for dispersion relation loaded')
k = Sdisp.k;
freq = Sdisp.freq;

%%
% points selected by hand 
kpts = [0.5443 0.5528];
omegapts = [2*pi*0.3878 2*pi*0.4097];
g = 9.81;
k_list = linspace(0.01,10,100); % array of k for plots
deep_water = sqrt(g*k_list);

figure, 
loglog(k,2*pi*freq,'o')
grid on 
hold on 
loglog(kpts,omegapts,'d')
hold on 
loglog(k_list,deep_water,'k-')
axis([0.06 2 0.8 6])
xlabel('$k \: \rm (m^{-1})$','Interpreter','latex')
ylabel('$\omega \: \rm (s^{-1})$','Interpreter','latex')
legend('','','$\omega =  \sqrt{gk}$','Location','southeast')
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;


% ################################
%% Get A(f,k) colorplot for Vx
% ################################
add_pow2 = [0 ,0 ,0];
facq_t = 1/m.ft;
facq_x = 1/m.fx;

% restriction to continuous ice
xmin = 1; %min index along x-direction
xmax = 48; %max index along x-direction
V = Vx(xmin:xmax,:,:);

[E,f,k,shift] = get_A_fk(V,add_pow2,facq_t,facq_x);

%%
save_name = [fig_folder 'A_fk_data_continuous_Vx_xmin_1xmax_48'];
save(save_name,'E','f','k','shift','-v7.3')
disp('Data Saved')
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
ylabel('$\omega \: \rm (Hz)$')
colormap(jet)
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'ColorScale','log')
cbar = colorbar();
cbar.Label.String = '$|\hat{V}_x|(k,\omega)$';
cbar.Label.Interpreter = 'latex';
hold on 
loglog(k_list,deep_water,'w--')
axis([0.1 2, 0.9 6])
set_Papermode(gcf)
set(gca,'FontSize',13)
caxis([0.9e-5 7e-3])
% hold on 
% loglog(k_list,yth,'w--')
% hold on
% loglog(k_list,harmonic2,'w--')
% hold on 
% loglog(k_list,harmonic3,'w--')

% lgnd = legend('',['$\omega^2 = gk \tanh(gh_w) \: h_w = ' num2str(h_w) '\: \rm m$'],'','');
% set(lgnd,'Location','southeast')
% set(lgnd,'color','none')
% loglog(k_list,harmonic2,'w--')
% hold on 
% loglog(k_list,harmonic3,'w--')
% colormap(slanCM('thermal-2'))

%% Plot FFT space for a given frequency
kx = -2*pi*facq_x*(-size(shift,1)/2:size(shift,1)/2-1)/size(shift,1);
ky = 2*pi*facq_x*(-size(shift,2)/2:size(shift,2)/2-1)/size(shift,2);

selected_freq = 0.42;
[~,idx] = min(abs(f - selected_freq));

A = abs(shift(:,:,idx));
figure, 
pcolor(kx,ky,A')
% shading interp
cbar = colorbar();
% set(gca,'ColorScale','log')
% caxis([1e-5 1e-2])

% ##############################
%% Work on velocity projection 
% ##############################
%% Select a frequency and the associated demodulated field 
% we also compute the angle of waves propagation : theta, for Vx and Vy 

% select a frequency
selected_freq = 0.42; % selected frequencies between which we proceed to the analysis
[min_freq, i0] = min(abs(f - selected_freq));

disp(['Selected frequency : ' num2str(f(i0)) ' Hz'])

%% Select the area to study (continuous or fragmented)
figure, 
% for continuous ice
xmin = 1; %in pixels
[ ~,xmax] = min(abs(36 - m.x)); %in pixels

% whole field 
% xmin = 1;
% xmax = size(m.Vx,1);

pcolor(m.x(xmin:xmax),m.y,real(FFT_Vx(xmin:xmax,:,i0))')
axis image
shading interp
colormap(redblue)

%%

padding_bool = 1;
add_pow2 = 2; % additional power of 2 for padding 
black_mask = 10;
caxis_amp = 0.002;

% show the demodulated field 
A_Vx = FFT_Vx(xmin:xmax,:,i0);
figure, 
pcolor(x_scaled(xmin:xmax),y_scaled,real(A_Vx)')
shading interp
axis image 
colormap(redblue)
colorbar()
caxis([-caxis_amp caxis_amp])

[~,kx_peak,ky_peak] = get_single_wavevector(A_Vx,padding_bool,add_pow2,fx,black_mask,-1);
% convert (kx,ky) -> (k,theta)
[theta_Vx,k_Vx] = cart2pol(kx_peak,ky_peak); 

disp(['Theta = ' num2str(theta_Vx) ' rad'])
disp(['k = ' num2str(k_Vx) ' m-1'])
A_Vy = FFT_Vy(xmin:xmax,:,i0);
figure, 
pcolor(x_scaled(xmin:xmax),y_scaled,real(A_Vy)')
shading interp
axis image 
colormap(redblue)
colorbar()
caxis([-caxis_amp caxis_amp])

[~,kx_peak,ky_peak] = get_single_wavevector(A_Vy,padding_bool,add_pow2,fx,black_mask,-1);
% convert (kx,ky) -> (k,theta)
[theta_Vy,k_Vy] = cart2pol(kx_peak,ky_peak); 

disp(['Theta = ' num2str(theta_Vy) ' rad'])
disp(['k = ' num2str(k_Vy) ' m-1'])
% final value of theta 
theta = mean([theta_Vx,theta_Vy]);

%% Projection of the velocity field : 

V_par = A_Vx*cos(theta) - A_Vy*sin(theta);
V_perp = A_Vx*sin(theta) + A_Vy*cos(theta);

% V_par = -A_Vx*cos(theta) - A_Vy*sin(theta);
% V_perp = A_Vx*sin(theta) - A_Vy*cos(theta);

figure, 
pcolor(m.x(xmin:xmax),m.y,real(V_par)')
shading interp
axis image 
colormap(redblue)
colorbar()
caxis([-caxis_amp caxis_amp])
xlabel('$x \: \rm (m)$','Interpreter','latex')
ylabel('$y \: \rm (m)$','Interpreter','latex')

figure,
pcolor(m.x(xmin:xmax),m.y,real(V_perp)')
shading interp
axis image 
colormap(redblue)
colorbar()
caxis([-caxis_amp caxis_amp])
xlabel('$x \: \rm (m)$','Interpreter','latex')
ylabel('$y \: \rm (m)$','Interpreter','latex')

%% 2D FFT on the parallel velocity

[~,kx_peak,ky_peak] = get_single_wavevector(V_par,padding_bool,add_pow2,fx,black_mask,-1);
[theta_Vpar,k_Vpar] = cart2pol(kx_peak,ky_peak); 
disp(['k = ' num2str(k_Vpar) ' m-1'])

%% Video of projected velocity field for all frequencies 

fx = 1/m.fx; % scale in pix / meter
selected_freq = [0.1 0.7];
x_bound = [1 48]; % corresponds to ice 
caxis_amp = -1;
left_bool = 1;
fig_name = 'Parallel_Demodulated_field_ice';
save_image = 0;
save_video = 0;

[~,idx_start] = min(abs(f - selected_freq(1)));
[~,idx_end] = min(abs(f - selected_freq(2)));
xmin = x_bound(1);
xmax = x_bound(2);

relevant_indices = (idx_start:1:idx_end);
clear T

fig_FFT_map = figure;
fig_FFT_map.Color = [1, 1, 1,];
for i = 1:length(relevant_indices)
    i0 = relevant_indices(i);
    A_Vx = FFT_Vx(xmin:xmax,:,i0);
    [~,kx_peak,ky_peak] = get_single_wavevector(A_Vx,padding_bool,add_pow2,fx,black_mask,-1);
    % convert (kx,ky) -> (k,theta)
    [theta_Vx,k_Vx] = cart2pol(kx_peak,ky_peak); 

    disp(['Theta = ' num2str(theta_Vx) ' rad'])
    disp(['k = ' num2str(k_Vx) ' m-1'])
    A_Vy = FFT_Vy(xmin:xmax,:,i0);

    [~,kx_peak,ky_peak] = get_single_wavevector(A_Vy,padding_bool,add_pow2,fx,black_mask,-1);
    % convert (kx,ky) -> (k,theta)
    [theta_Vy,k_Vy] = cart2pol(kx_peak,ky_peak); 

    disp(['Theta = ' num2str(theta_Vy) ' rad'])
    disp(['k = ' num2str(k_Vy) ' m-1'])
    % final value of theta 
    theta = mean([theta_Vx,theta_Vy]);

    % Projection of the velocity field : 
    V_par = A_Vx*cos(theta) - A_Vy*sin(theta);
    V_perp = A_Vx*sin(theta) + A_Vy*cos(theta);
        
    if save_video
        video_filename = [fig_folder fig_name '.avi']; % folder where the video is saved
        vid = VideoWriter(video_filename);
        vid.FrameRate = 3;
        open(vid)
    end 

    pcolor(m.x(xmin:xmax),m.y,real(V_par)')
    shading interp
    axis image 
    colormap(redblue)
    colorbar()
    xlabel('$x \: \rm (m)$','Interpreter','latex')
    ylabel('$y \: \rm (m)$','Interpreter','latex')
    
    if caxis_amp > 0 
        caxis([-caxis_amp caxis_amp])
    end 
    ax = gca;
    ax.FontSize = 13;
    
    getframe();
    frequency = f(i0);
    if save_video
        title(['$f = ' num2str(frequency) ' \: \rm(Hz)$'],'Interpreter','latex')
    end 
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'filename','-dpdf','-r0')
    %pause(0.2)
    
    if save_image
        frequency_p = strrep(num2str(frequency),'.','p'); % replace point by p
        fig_file_fig = [new_folder 'FFT_map_f_' frequency_p];
        saveas(fig_FFT_map,fig_file_fig,'fig')
    end
    
    if save_video
        T(i) = getframe(gcf);     
    end 
    
end

% Write video
if save_video
    all_valid = true;
    flen = length(T);
    for K = 1 : flen
      if isempty(T(K).cdata)
        all_valid = false;
        fprintf('Empty frame occurred at frame #%d of %d\n', K, flen);
      end
    end
    if ~all_valid
       error('Did not write movie because of empty frames')
    end

    writeVideo(vid,T)
    close(vid)
end 

%% Projection of velocity field + 2D FFT 

fx = 1/m.fx; % scale in pix / meter
selected_freq = [0.1 0.7];
x_bound = [1 48]; % corresponds to ice 
caxis_amp = -1;
padding_bool = 1;
add_pow2 = 2;
black_mask = 10;
fig_name = 'Parallel_Demodulated_field_ice';
save_image = 0;
save_video = 0;

[~,idx_start] = min(abs(f - selected_freq(1)));
[~,idx_end] = min(abs(f - selected_freq(2)));
xmin = x_bound(1);
xmax = x_bound(2);

relevant_indices = (idx_start:1:idx_end);
V_par = zeros(size(FFT_Vx(xmin:xmax,:,idx_start:idx_end)));
V_perp = zeros(size(FFT_Vx(xmin:xmax,:,idx_start:idx_end)));
wave_vector = zeros(1,length(relevant_indices));
freq = f(relevant_indices);

for i = 1:length(relevant_indices)
    i0 = relevant_indices(i);
    A_Vx = FFT_Vx(xmin:xmax,:,i0);
    [~,kx_peak,ky_peak] = get_single_wavevector(A_Vx,padding_bool,add_pow2,fx,black_mask,-1);
    % convert (kx,ky) -> (k,theta)
    [theta_Vx,k_Vx] = cart2pol(kx_peak,ky_peak); 

    disp(['Theta = ' num2str(theta_Vx) ' rad'])
    disp(['k = ' num2str(k_Vx) ' m-1'])
    A_Vy = FFT_Vy(xmin:xmax,:,i0);

    [~,kx_peak,ky_peak] = get_single_wavevector(A_Vy,padding_bool,add_pow2,fx,black_mask,-1);
    % convert (kx,ky) -> (k,theta)
    [theta_Vy,k_Vy] = cart2pol(kx_peak,ky_peak); 

    disp(['Theta = ' num2str(theta_Vy) ' rad'])
    disp(['k = ' num2str(k_Vy) ' m-1'])
    % final value of theta 
    theta = mean([theta_Vx,theta_Vy]);

    % Projection of the velocity field : 
    V_par(:,:,i) = A_Vx*cos(theta) - A_Vy*sin(theta);
    V_perp(:,:,i) = A_Vx*sin(theta) + A_Vy*cos(theta);
    
    % 2D FFT of the parallel velocity field 
    [~,kx_peak,ky_peak] = get_single_wavevector(V_par(:,:,i),padding_bool,add_pow2,fx,black_mask,-1);
    [~,k] = cart2pol(kx_peak,ky_peak);
    wave_vector(i) = k;
    
end 

%%
figure, 
loglog(wave_vector,freq,'o')
grid on 


%% Try to apply a hamming function
img = A_Vx;

[Nx, Ny] = size(img);
Nb_elements = Nx*Ny; % number of elements of matrix a

if padding_bool 
        % generate a 2D - hamming function
        ham_window = window2(Nx,Ny,@hamming);
        
        padding_x = 2^(nextpow2(Nx) + add_pow2);
        padding_y = 2^(nextpow2(Ny) + add_pow2);
        fft_2D = fft2(img.*ham_window,padding_x,padding_y);
        disp('Padding used')
else 
        fft_2D = fft2(img);
        padding_x = Nx;
        padding_y = Ny;
        disp('No padding')
end

fft_2D = fft_2D/Nb_elements; % normalization by the initial number of elements
% #################################
% SCALING : create kx and ky arrays
% #################################
kx = -2*pi*fx*(-padding_x/2:padding_x/2-1)/padding_x;
ky = 2*pi*fx*(-padding_y/2:padding_y/2-1)/padding_y;

shifted_fft = fftshift(fft_2D);

%%

figure, 
pcolor(kx,ky,abs(shifted_fft)')
shading interp
axis image
axis([-1.5 1.5 -1.5 1.5])
xlabel('$k_x$','Interpreter','latex');
ylabel('$k_y$','Interpreter','latex');
cbar = colorbar();
cbar.Label.String = '$ | \hat{V_x}| (k_x,k_y,f)$';

ax = gca;
ax.FontSize = 13;
%% Do a movie of the velocity field 

save_video = 0;

if save_video
    video_filename = [fig_folder 'Velocity_field_Vy_video.avi']; % folder where the video is saved
    vid = VideoWriter(video_filename);
    vid.FrameRate = 30;
    open(vid)
end 

[nx,ny,nt] = size(m.Vx);
x = (1:1:nx);
y = (ny:-1:1);
[X,Y] = meshgrid(x,y);

velocity_fig = figure;
velocity_fig.Color = [1 , 1 ,1];

for i = 1 : nt
%     V = sqrt(m.Vx(:,:,i).^2 + m.Vy(:,:,i).^2);
    
    pcolor(x/fx,y/fx,m.Vy(:,:,i)')
    shading interp
    title(['Frame ' num2str(i)])
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    ax = gca;
    ax.FontSize = 13;
    axis image
    colorbar()
    caxis([-3 3])
    if save_video
        T(i)=getframe(gcf);     
    end
    pause(0.1)

end 

if save_video
    all_valid = true;
    flen = length(T);
    for K = 1 : flen
      if isempty(T(K).cdata)
        all_valid = false;
        fprintf('Empty frame occurred at frame #%d of %d\n', K, flen);
      end
    end
    if ~all_valid
       error('Did not write movie because of empty frames')
    end

    writeVideo(vid,T)
    close(vid)    
end 


%% ############ FUNCTION SECTION ##########

function [A_peak,kx_peak,ky_peak] = get_single_wavevector(A,padding_bool,add_pow2,fx,black_mask,caxis_amp)
% this function computes the 2D-FFT of a 2D field, apply a mask at the
% center of the 2D-FFT, and select the associated coordinates of the wave
% vector 

% It takes as arguments : 
% - A : the 2D field [nx,ny], it can be complex
% - padding_bool : a boolean to choose to perform padding
% - add_pow2 : int, if additional padding is needed 
% - fx : the acquisition frequency (in pixels/meters)
% - black_mask : interger that defines the size of the window, centered on
% zero, to be set to zero.
% - caxis_amp : value to set the colobar axis 

    % computes the 2D FFT
    [shifted_fft,fft_2D,kx,ky] = spatial_FFT(A,padding_bool,add_pow2,fx);
    % kx = dimension 1
    % ky = dimension 2

    % Mask center 
    if black_mask >= 0 % if black_mask < 0, we don't mask the center of the 2D-FFT
        shifted_fft(size(shifted_fft,1)/2+1-black_mask:size(shifted_fft,1)/2+1+black_mask,size(shifted_fft,2)/2+1-black_mask:size(shifted_fft,2)/2+1+black_mask)=zeros;
    end

    % Try to mask bottom 
    %shifted_fft(size(shifted_fft,1)/2 + 30:end,:) = zeros;
%     figure,
%     pcolor(kx,ky,abs(shifted_fft)');
%     shading interp
%     %imagesc(abs(shifted_fft));
% 
%     xlabel('$k_x$','Interpreter','latex');
%     ylabel('$k_y$','Interpreter','latex');
%     cbar = colorbar();
%     cbar.Label.String = '$ | \hat{V_x}| (k_x,k_y,f)$';
%     cbar.Label.Interpreter = 'latex';
%     if caxis_amp > 0
%         caxis([0 caxis_amp])
%     end 
%     %axis([-4 4 -4 4])
%     axis([-1.5 1.5 -1.5 1.5]);
%     ax = gca;
%     ax.FontSize = 13;

    [A_peak, peak_position] = max(abs(shifted_fft), [],'all', 'linear');
    % convert linear index to subscript
    sz = size(shifted_fft); % size of the matrix 2D FFT
    [row, col] = ind2sub(sz,peak_position); % position of the peak in the matrix
    hold on 
    % plot position of the peak
    kx_peak = kx(row); % dimension 1
    ky_peak = ky(col); % dimension 2

    plot(kx_peak,ky_peak,'ro');
    hold off

end 
