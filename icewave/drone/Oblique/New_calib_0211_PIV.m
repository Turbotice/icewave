clear all; 
close all; 

%% Load PIV 

base = 'E:/Rimouski_2024/Data/2024/0211/Drones/bernache/matData/18-stereo_001/';
%base = 'E:/PIVlab_drone/matdata/DJI_0308_Dt4_W64_full/';

filename = 'PIV_processed_i011496_N15496_Dt4_b1_W32_full_total_processed.mat';
%%
disp('Loading data..')
load([base filename])
disp('Data loaded')

%% Create an inclination parameter
param.alpha_0 = m.alpha_0;
param.H = m.h_drone; 
param.unit_alpha_0 = 'deg';
param.unit_H = 'meter';

param_filename = [base 'Param_Debut_mesange_PIV_DJI_20240211213233_0222_D.mat'];
save(param_filename,'param','-v7.3')

%% Set folder where to save plots 

fig_folder =  [base 'inclination_plots/'];

if ~exist(fig_folder)
    mkdir(fig_folder)
end


%%

H = param.H; % Hauteur en mètres
alpha_0 = param.alpha_0*pi/180; %angle passé en radians
theta_x = 34.15;

L_x = 3840; % size of the image in pixel, along x-axis
fx_pix = L_x/(2*H*tan(theta_x*pi/180)); % scale in pixels / meter
W = 32; 

f = 2830; % focale apparente en pixels 
%%
Vy = m.Vy;
Vx = m.Vx;
Dt = 4;
fps = 29.97;

dt = Dt/fps;
x_pix = repmat([1:1:size(Vy,1)]*W/2 + 1,size(Vy,2),1)';
y_pix = repmat([1:1:size(Vy,2)]'*W/2 + 1,1,size(Vy,1))';

%% 
% Les distances sont à prendre par rapport au centre du capteur

x_0 = 3840/2;
y_0 = 2160/2;

%% Definition of X and Y in real framework

Y = (y_pix - y_0)*H/sin(alpha_0)./(f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0));
X = (x_pix - x_0)*H./(f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0));

X = X;
Y = -Y;

%% Check if redressement is correct

figure, 
for i = 1:size(Vy,3)
    
    surf(X,Y,Vy(:,:,i))
    shading interp
    view(2)
    title(num2str(i))
    pause(0.01)
    colormap(redblue)
end


%% Computation of vertical displacement, test #3

vid_name = 'Vertical_displacement_mesange';

save_video = 0;

if save_video
    video_filename = [fig_folder vid_name '.avi']; % folder where the video is saved
    vid = VideoWriter(video_filename);
    vid.FrameRate = fps;
    open(vid)
end


for i0 = 1:1:size(Vy,3)
    
    vx = Vx(:,:,i0);
    vy = Vy(:,:,i0);
    
%     y_prime = y_pix + vy;
%     x_prime = x_pix + vx;
    
    Vz = H*f*vy./((y_pix - y_0)*cos(alpha_0) + f*sin(alpha_0))...
    ./((y_pix - y_0 + vy)*sin(alpha_0) - f*cos(alpha_0));

    dz(:,:,i0) = Vz.*fps/Dt;
    
    disp(i0)
    
    surf(X,Y,squeeze(dz(:,:,i0)))
    shading interp
    cbar = colorbar();
    cbar.Label.String = '$V_z \: \rm (m/s)$';
    cbar.Label.Interpreter = 'latex';
    colormap(redblue)
    caxis([-2.5 2.5])
    view(2)
    title(num2str(i0))
    
    if save_video
        T(i0) = getframe(gcf);
    else 
        pause(0.01)
    end 
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

%% Positions bouées + téléphones
%B4 : 2594 x 750 % bernache  %B4 : 1126 x 891 % mesange
%B5 : 3204 x 759 % bernache  %B5 : 586 x 911 % mesange

x_b = 586; % x-position in pixel 
y_b = 911; % y-position in pixel

% X and Y coordinate in sea ice frame work
Y_b = (y_b - y_0)*H/sin(alpha_0)./(f*sin(alpha_0) + (y_b - y_0).*cos(alpha_0));
X_b = (x_b - x_0)*H./(f*sin(alpha_0) + (y_b - y_0).*cos(alpha_0));

Y_b = - Y_b;

% Show position of the studied point
map_fig = figure; 
plot3(X_b,Y_b,100,'or');
hold on 
surf(X,Y,Vy(:,:,1))
shading interp
view(2)
title('Buoy #5')

figname = [fig_folder 'Map_Buoy_5_minus_Vy'];
saveas(map_fig,figname,'fig')
saveas(map_fig,figname,'pdf')

%%
% find closest box to the real position of the POI
idx_x = floor(x_b*2/W);
idx_y = floor(y_b*2/W);

pixel_fig = figure; 
plot3(X(idx_x,idx_y),Y(idx_x,idx_y),100,'db');
hold on 
surf(X,Y,Vy(:,:,1))
hold on 
plot3(X_b,Y_b,100,'or');
shading interp
view(2)
title('Map buoy #5')

figname = [fig_folder 'Closest_pixel_Buoy_5_minus_Vy'];
saveas(pixel_fig,figname,'fig')
saveas(pixel_fig,figname,'pdf')

%%

% img_path = ;
% img = imread(img_path);
% 
% figure, 
% surf(X,Y,img)
% hold on 

x_c = floor(size(dz,2)/2);
y_c = floor(size(dz,2)/2 + 1);
% x_c = 90;
% y_c = 40;
raw_fig = figure; 
plot(m.t,squeeze(-dz(idx_x,idx_y,:)))
% hold on 
% plot(m.t,squeeze(dz_2(x_c,y_c,:)))
xlabel('$t \: \rm (s)$','Interpreter','latex')
ylabel('$V_z \: \rm (m/s)$','Interpreter','latex')
% legend('Method Antonin','Method beta-angle')
grid on

figname = [fig_folder 'Raw_time_plot_buoy_5_minus_Vy'];
saveas(raw_fig,figname,'fig')
saveas(raw_fig,figname,'pdf')

%% Apply a low pass filter to the data

signal = squeeze(-dz(idx_x,idx_y,:));
fpass = 2; % cutoff frequency of the low pass filter
passed = lowpass(signal,fpass,fps); % apply the low pass filter

low_passed_fig = figure; 
plot(m.t,passed)
xlabel('$t \: \rm (s)$','Interpreter','latex')
ylabel('$V_z \: \rm (m/s)$','Interpreter','latex')
grid on 

figname = [fig_folder 'Low_pass_fcut_' num2str(fpass) '_time_plot_buoy_5_minus_Vy'];
saveas(low_passed_fig,figname,'fig')
saveas(low_passed_fig,figname,'pdf')

%% Do a spatio-temporal 

signal2 = squeeze(-dz(:,idx_y,:));

spatio_temp_fig = figure; 
imagesc(X(:,idx_y),m.t,signal2')
colorbar()
ylabel('$t \: \rm (s)$','Interpreter','latex')
xlabel('$x \: \rm (m)$','Interpreter','latex')

figname = [fig_folder 'Saptiotemporal_buoys_line_minus_Vy'];
saveas(spatio_temp_fig,figname,'fig')
saveas(spatio_temp_fig,figname,'pdf')

%% Create a structure for reconstruction

% save main fields
s.Vx = Vx;
s.Vy = Vy;
s.dz = -dz;
s.X = X;
s.Y = Y;
s.t = m.t;
s.drone = 'mesange';
% save main parameters
s.param = struct('H',H,'alpha_0',alpha_0,'focale',f,'Dt',Dt,'fps',fps,'W',W);
s.unit_param = struct('unit_H','meter','unit_alpha_0','rad','unit_f','pix');

%%

matname = [base 'Data_PIV_oblique_Dt' num2str(Dt) '_W' num2str(W) '_mesange_minus_Vy'];
save(matname,'s','-v7.3')

% ############################################
%% Load reconstructed data from both drones
% ############################################
base_bernache = 'W:/SagWin2024/Data/0211/Drones/bernache/matData/18-stereo_001/';
base_mesange = 'W:/SagWin2024/Data/0211/Drones/mesange/matData/2-stereo_001/';

file_bernache = [base_bernache 'Data_PIV_oblique_Dt4_W32_bernache.mat'];
file_mesange = [base_mesange 'Data_PIV_oblique_Dt4_W32_mesange.mat'];

disp('Loading data...')
s_bernache = load(file_bernache);
% s_mesange = load(file_mesange);
disp('Data loaded')

%% Create a MP4 movie
fig_folder = 'W:/SagWin2024/Data/0211/Drones/bernache/matData/18-stereo_001/';
vid_name = 'Vertical_displacement_bernache_redblue_caxis_1p5';

save_video = 1;
X = s_bernache.s.X;
Y = s_bernache.s.Y;
dz = s_bernache.s.dz;
fps = s_bernache.s.param.fps;

if save_video
    video_filename = [fig_folder vid_name]; % folder where the video is saved
    vid = VideoWriter(video_filename,'MPEG-4');
    vid.FrameRate = fps;
    vid.Quality = 90;
    open(vid)
end


for i0 = 1:1:size(dz,3)
    
    surf(X,Y,squeeze(dz(:,:,i0)))
    shading interp
    cbar = colorbar();
    cbar.Label.String = '$V_z \: \rm (m/s)$';
    cbar.Label.Interpreter = 'latex';
    xlabel('$x \: \rm (m)$')
    ylabel('$y \: \rm (m)$')
    colormap(redblue)
    caxis([-1.5 1.5])
    view(2)
    title(['$t = ' num2str(s_bernache.s.t(i0),'%4.2f') ' \: \rm s$'])

    ax = gca;
    ax.FontSize = 13;
    set_Papermode(gcf)
    if save_video
        T(i0) = getframe(gcf);
    else 
        pause(0.01)
    end 
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
%% Plot both signal for buoy #4
fig_folder = 'E:/Rimouski_2024/Data/2024/0211/Drones/bernache_mesange_superposition_plots/' ;

if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Plot both signal for a given object
%# Buoy #5
obj_txt = 'buoy_5';
x_bernache = 3204;
y_bernache = 759;

x_mesange = 586;
y_mesange = 911;

%# Buoy #4
% obj_txt = 'buoy_4';
% x_bernache = 2594;
% y_bernache = 750;
% 
% x_mesange = 1126;
% y_mesange = 891;

%# Buoy #3 
% obj_txt = 'buoy_3';
% x_bernache = 2030;
% y_bernache = 756;
% 
% x_mesange = 1812;
% y_mesange = 857;

signal_bernache = get_1D_signal(x_bernache,y_bernache,s_bernache.s.dz,W);
signal_mesange = get_1D_signal(x_mesange,y_mesange,s_mesange.s.dz,W);
f_cut = 2; % cut off frequency, 2 Hz
low_bernache = lowpass(signal_bernache,f_cut,fps);
low_mesange = lowpass(signal_mesange,f_cut,fps);

time_array = s_bernache.s.t;
idx_min = 1;
idx_max = length(time_array);

filter_fig = figure;
plot(time_array(idx_min:idx_max),low_bernache(idx_min:idx_max))
hold on 
plot(time_array(idx_min:idx_max),low_mesange(idx_min:idx_max))
grid on 
legend('bernache','mesange')
xlabel('$t \: \rm (s)$', 'Interpreter','latex')
ylabel('$V_z(t) \: \rm (m.s^{-1})$','Interpreter','latex')

raw_fig_superposition = figure; 
plot(time_array(idx_min:idx_max),signal_bernache(idx_min:idx_max))
hold on 
plot(time_array(idx_min:idx_max),signal_mesange(idx_min:idx_max))
grid on 
legend('bernache','mesange')
xlabel('$t \: \rm (s)$', 'Interpreter','latex')
ylabel('$V_z(t) \: \rm (m.s^{-1})$','Interpreter','latex')

save_bool = 0;

if save_bool

    figname = [fig_folder 'Supperposition_' obj_txt '_raw_datas'];
    saveas(raw_fig_superposition,figname,'fig')
    saveas(raw_fig_superposition,figname,'pdf')

    figname = [fig_folder 'Supperposition_' obj_txt '_filtered_datas'];
    saveas(filter_fig,figname,'fig')
    saveas(filter_fig,figname,'pdf')
end
%% Load raw velocity fields 

base_bernache = 'E:/Rimouski_2024/Data/2024/0211/Drones/bernache/matData/18-stereo_001/';
base_mesange = 'E:/Rimouski_2024/Data/2024/0211/Drones/mesange/matData/2-stereo_001/';

file_bernache = [base_bernache 'PIV_processed_i011500_N15500_Dt4_b1_W32_full_total_processed.mat'];
file_mesange = [base_mesange 'PIV_processed_i011496_N15496_Dt4_b1_W32_full_total_processed.mat'];

disp('Loading data...')
raw_bernache = load(file_bernache);
raw_mesange = load(file_mesange);
disp('Data loaded')

%% Norm of velocity fields

V_bernache = sqrt(raw_bernache.m.Vx.^2 + raw_bernache.m.Vy.^2);
V_mesange = sqrt(raw_mesange.m.Vx.^2 + raw_mesange.m.Vy.^2);

% ###############################################################
%% ################## Load Buoys signal #########################
% ###############################################################

base_buoys = 'W:/SagWin2024/Data/0211/BoueeVague/B5/mat/';

buoy_filename = [base_buoys 'buoy5_sbg_20240211_2000.mat'];

disp('Loading data..')
load(buoy_filename)
disp('Data loaded')

%% Plot raw data from buoys
buoy_vz = IMU.IMU_DATA.ACCEL_Z ; % get buoy acceleration
buoy_UTC = IMU.IMU_DATA.t ;
%buoy_UTC = datetime(buoy_UTC, 'convertfrom', 'posixtime', 'Format', 'MM/dd/yy HH:mm:ss.SSS');

figure, plot(buoy_UTC,buoy_vz)

% integrate to get velocity 

%% Select the UTC-time associated to the beginning of the acquisition

Y = IMU.UTC_TIME.YEAR(1);
M = IMU.UTC_TIME.MONTH(1);
D = IMU.UTC_TIME.DAY(1);
H = IMU.UTC_TIME.HOUR(1);
MIN = IMU.UTC_TIME.MIN(1);
S = IMU.UTC_TIME.SEC(1);
% MS = IMU.UTC_TIME.NANOSEC(1);

% initial time of recording for buoys
t0 = datetime(Y,M,D,H,MIN,S,'TimeZone','UTC'); 
t0.TimeZone = 'UTC'; % converts time to UTC time 
t0.Format = 'yyyy-MMM-dd HH:mm:ss';

%% Convert this UTC time to time since epoch 

t0_epoch = convertTo(t0,'posixtime');
% create an array of time t since the epoch for the acquisition
t_epoch = t0_epoch + IMU.IMU_DATA.t ; 

% create an array of UTC time for the acquisition 
UTC_t = datetime(t_epoch, 'convertfrom', 'posixtime', 'Format', 'MM/dd/yy HH:mm:ss.SSS');
UTC_t.TimeZone = 'UTC';
%% Plot buoy data with epoch time and UTC time 
figure, 
plot(t_epoch,IMU.IMU_DATA.ACCEL_Z)

figure, 
plot(UTC_t,IMU.IMU_DATA.ACCEL_Z,'.')
grid on 

%% Plot bernache and mesange signal with UTC_time 
%# Buoy #1
% obj_txt = 'buoy_1';
% x_bernache = 1284;
% y_bernache = 752; 
% 
% x_mesange = 2570;
% y_mesange = 821;

%# Buoy #2
% obj_tkt = 'buoy_2';
% x_bernache = 1457;
% y_bernache = 738;
% 
% x_mesange = 2401;
% y_mesange = 846;

%# Buoy #3 
% obj_txt = 'buoy_3';
% x_bernache = 2029;
% y_bernache = 755;
% 
% x_mesange = 1812;
% y_mesange = 857;

%# Buoy #5
obj_txt = 'buoy_5';
x_bernache = 3204;
y_bernache = 759;

x_mesange = 586;
y_mesange = 911;

%# Buoy #4
% obj_txt = 'buoy_4';
% x_bernache = 2593;
% y_bernache = 749;
% 
% x_mesange = 1126;
% y_mesange = 891;

signal_bernache = get_1D_signal(x_bernache,y_bernache,s_bernache.s.dz,32);
signal_mesange = get_1D_signal(x_mesange,y_mesange,s_mesange.s.dz,32);

figure, 
plot(s_bernache.s.UTC_time',signal_bernache)
hold on 
plot(s_bernache.s.UTC_time',signal_mesange)
grid on 

%% select relevant signal 
t_start = datetime(2024,2,11,20,34,45,'TimeZone','UTC'); 
t_end = datetime(2024,2,11,20,37,26,'TimeZone','UTC');

mask = (UTC_t >= t_start) & (UTC_t < t_end);
ROI_signal = IMU.IMU_DATA.ACCEL_Z(mask);
ROI_t = UTC_t(mask);

figure, 
plot(ROI_t,ROI_signal)

%% Applies a filter and Computes vertical velocity of boys
fs = 50; %acquisition frequency of buoys
fc = 0.1; %cutoff frequency 

[b,a] = butter(6,fc/(fs/2),"high");
% [b,a] = butter(6,/25,"high");
filtered_buoy = filtfilt(b,a,double(ROI_signal));

figure, 
plot(ROI_t,filtered_buoy)
grid on 

%% Select new ROI signal 

t_start = datetime(2024,2,11,20,34,50,'TimeZone','UTC'); 
t_end = datetime(2024,2,11,20,37,26,'TimeZone','UTC');

mask = (ROI_t >= t_start) & (ROI_t < t_end);
filtered_buoy = filtered_buoy(mask);
ROI_t = ROI_t(mask);

figure, 
plot(ROI_t,filtered_buoy)


%% Compute the vertical velocity

% dt = mean(diff(ROI_t)); %time step
vz = cumsum(filtered_buoy - mean(filtered_buoy))/fs;
figure, 
plot(ROI_t,vz)

%% Superposition of drones and buoy with UTC time 
delta_t0 = s_bernache.s.UTC_time(1) - ROI_t(1); % difference of t0 between bernache UTC time and buoys UTC time 
add_dt = milliseconds(500); % add 500ms 
UTC_drone = s_bernache.s.UTC_time - delta_t0 + add_dt;

figure, 
plot(UTC_drone',signal_bernache)
hold on 
plot(UTC_drone',signal_mesange)
hold on 
plot(ROI_t,vz,'k','LineWidth',1.0)
grid on 
xlabel('$t$','Interpreter','latex')
ylabel('$V_z \: \rm (m.s^{-1})$','Interpreter','latex')
title('Buoy 2','Interpreter','latex')
legend('Bernache','Mesange','Buoy','Interpreter','latex')
ylim([-1 1])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;

% ######################################
%% Get Drone signal noise 
% ######################################

%% Set fig_folder 
fig_folder =  ['W:/SagWin2024/Data/0211/Drones/', 'Drone_noise/'];

if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Create a matrix with all signals for a given buoy 
B(1) = struct('signal',vz,'UTC_t',ROI_t,'fps',vz);
B(2) = struct('signal',signal_bernache,'UTC_t',UTC_drone','fps',s_bernache.s.param.fps);
B(3) = struct('signal',signal_mesange,'UTC_t',UTC_drone','fps',s_bernache.s.param.fps);

%% Interpolate Buoy signal with a sampling similar to drone sampling 
I1 = interp1(B(1).UTC_t,B(1).signal,B(2).UTC_t);
% Compare interpolated signal with original one 
figure, 
plot(B(2).UTC_t,I1,'.')
hold on 
plot(B(1).UTC_t,B(1).signal)
grid on 

%% FFT of each signal 
padding_bool = 0;
add_pow2 = 0;
for i = 1:3
    if i > 1
        signal = B(i).signal;
        [TF,TF_pos,f] = fft_1D(signal,B(i).fps,padding_bool,add_pow2);
        B(i).TF_pos = TF_pos;
        B(i).TF = TF;
        B(i).f = f;
    else
        signal = I1;
        [TF,TF_pos,f] = fft_1D(signal,B(2).fps,padding_bool,add_pow2);
        B(i).TF_pos = TF_pos;
        B(i).TF = TF;
        B(i).f = f;
    end 
end 

%% Superpose FFT of each signal 

figure, 
c = [0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980];
for i = 1:3
    loglog(B(i).f,abs(B(i).TF_pos),'Color',c(i,:))
    hold on 
end
grid on 

xlabel('$f \: \rm (Hz)$')
ylabel('$\hat{s}(f) \: \rm (m.s^{-1})$')
axis([5e-3 30 1e-6  1])
ax = gca;
ax.FontSize = 13;
set_Papermode(gcf)
legend('Buoy','Bernache','Mesange')

figname = [fig_folder 'Superposition_FFT_Buoy5'];
saveas(gcf,figname,'fig')
saveas(gcf,figname,'pdf')

%% Divide FFT of bernache / mesange by FFT of buoy signal (only absolute value)

norm_bernache = abs(B(2).TF_pos)./abs(B(1).TF_pos);
norm_mesange = abs(B(3).TF_pos)./abs(B(1).TF_pos);

B(2).TF_noise = B(2).TF./abs(B(1).TF);
B(3).TF_noise = B(3).TF./abs(B(1).TF);

figure, 
loglog(B(2).f,norm_bernache)
hold on 
loglog(B(2).f,norm_mesange)
grid on 

xlabel('$f \: \rm (Hz)$')
ylabel('$\hat{s}/\hat{s}_{buoy} (f) \: \rm (m.s^{-1})$')
axis([5e-3 30 1e-1  1e3])
ax = gca;
ax.FontSize = 13;
set_Papermode(gcf)
legend('Bernache','Mesange','Location','southwest')

figname = [fig_folder 'FFT_drone_noise_buoy5'];
saveas(gcf,figname,'fig')
saveas(gcf,figname,'pdf')

%% Inverse Fourier transform 
noise_bernache = ifft(B(2).TF_noise*numel(B(2).TF_noise));
noise_mesange = ifft(B(3).TF_noise*numel(B(3).TF_noise));

figure, 
plot(B(2).UTC_t,noise_bernache)
grid on 
xlabel('$\rm Time \: (UTC)$')
ylabel('$s \: \rm (m.s^{-1})$')
% axis([5e-3 30 1e-1  1e3])
ax = gca;
ax.FontSize = 13;
set_Papermode(gcf)

figname = [fig_folder 'Bernache_noise_buoy5'];
saveas(gcf,figname,'fig')
saveas(gcf,figname,'pdf')

%% Save Buoy Matrix 
matname = [fig_folder 'Signal_initial_buoys_position'];
save(matname,'B_matrix','-v7.3')

%% Fourier spectrum for several pixels of displacement field 

% select a buoy #1 -> buoy 1, #2 -> buoy 4, #3 -> buoy 5
% buoy_idx = 1;
%# Buoy #1
% obj_txt = 'buoy_1';
x_bernache = 1284;
y_bernache = 752; 

x_mesange = 2570;
y_mesange = 821;

% get closest boxes 
% # (1,:) -> i,j box on bernache 
% # (2,:) -> i,j box on mesange 
W = 32;
idx_buoy(1,:) = [floor(x_bernache*2/W) floor(y_bernache*2/W)];
idx_buoy(2,:) = [floor(x_mesange*2/W) floor(y_mesange*2/W)];

window = 3; % semi-width of window around the buoy
Vz_local(1,:,:,:) = s_bernache.s.dz(idx_buoy(1,1)-window : idx_buoy,:);

% ###################################
%% Plot data at different positions 
% ###################################
fig_folder = 'W:/SagWin2024/Data/0211/Drones/bernache_mesange_superposition_plots/';

%# Buoy #5
% obj_txt = 'buoy_5';
% x_bernache = 3204;
% y_bernache = 759;
% 
% x_mesange = 586;
% y_mesange = 911;

%# Buoy #4
% obj_txt = 'buoy_4';
% x_bernache = 2593;
% y_bernache = 749;
% 
% x_mesange = 1126;
% y_mesange = 891;

%# Buoy #3 
% obj_txt = 'buoy_3';
% x_bernache = 2029;
% y_bernache = 755;
% 
% x_mesange = 1812;
% y_mesange = 857;

%# Buoy #2
obj_tkt = 'buoy_2';
x_bernache = 1457;
y_bernache = 738;

x_mesange = 2401;
y_mesange = 846;

%# Buoy #1
% obj_txt = 'buoy_1';
% x_bernache = 1284;
% y_bernache = 752;
% 
% x_mesange = 2570;
% y_mesange = 821;

%# Tel #1
% obj_txt = 'tel_1';
% x_bernache = 1438;
% y_bernache = 1248;
% 
% x_mesange = 2299;
% y_mesange = 410;

signal_bernache = get_1D_signal(x_bernache,y_bernache,s_bernache.s.dz,32);
signal_mesange = get_1D_signal(x_mesange,y_mesange,s_mesange.s.dz,32);

time_fig = figure; 
plot(s_bernache.s.t,signal_bernache)
hold on 
plot(s_mesange.s.t,signal_mesange);
grid on 
legend('bernache','mesange')
xlabel('$t \: \rm (s)$','Interpreter','latex')
ylabel('$V_z \: \rm (m/s)$', 'Interpreter','latex')

set_Papermode(time_fig)
subfolder = [fig_folder 'Vz/' obj_txt '/'];

if ~exist(subfolder)
    mkdir(subfolder)
end

figname = ['Time_evolution_Vz_' obj_txt];
saveas(time_fig, [subfolder figname],'fig')
saveas(time_fig, [subfolder figname],'pdf')

%% Plot velocity maps for different frames index

i0 = 1000;
img_bernache = raw_bernache.m.Vy(:,:,i0);
img_mesange = raw_mesange.m.Vy(:,:,i0);

figure, 
imagesc(img_bernache')
axis image

figure, 
imagesc(img_mesange')
axis image

%% Create UTC-time array 
% t0 Mesange
Y = 2024;
M = 02;
D = 11;
H = 21;
MIN = 35;
S = 10;
MS = 589;

% initial time of recording
t0_mesange = datetime(Y,M,D,H,MIN,S,MS,'TimeZone','Europe/Paris'); 
t0_mesange.TimeZone = 'UTC'; % converts time to UTC time 
t0_mesange.Format = 'yyyy-MMM-dd HH:mm:ss.SSS';

% t0 Bernache
Y = 2024;
M = 02;
D = 11;
H = 15;
MIN = 35;
S = 11;
MS = 690;

% initial time of recording
t0_bernache = datetime(Y,M,D,H,MIN,S,MS,'TimeZone','America/Montreal'); 
t0_bernache.TimeZone = 'UTC'; % converts time to UTC time 
t0_bernache.Format = 'yyyy-MMM-dd HH:mm:ss.SSS';

%%
% converts time array (in second) to UTC time 

% convert time array to duration in seconds 
time_s = milliseconds(s_mesange.snew.t*1000); 
t0 = t0_mesange;
t0.Format = 'HH:mm:ss.SSS'; % displayed format
UTC_time_array = t0 + time_s; % array of UTC time 

%%
Dt = 4;
W = 32;
s = s_mesange.snew;
s.UTC_time = UTC_time_array;
%%
matname = [base_mesange 'Data_PIV_oblique_Dt' num2str(Dt) '_W' num2str(W) '_mesange'];
save(matname,'s','-v7.3')

%% FUNCTION SECTION

function [signal] = get_1D_signal(x_obj,y_obj,dz,W)

% This function computes the vertical velocity at a given position of the
% PIV field, and returns it as a 1D signal
% The function takes the following arguments : 
% - x_obj : x-coordinate on the camera (in pixels)
% - y_obj : y-coordinate on the camera (in pixels)
% - dz : vertical velocity matrix [nx,ny,nt]
% - W : window size used for the PIV processing (pixels)

% x_0 = 3840/2;
% y_0 = 2160/2;

% object position in the sea ice framework
% Y_b = (y_obj - y_0)*H/sin(alpha_0)./(f*sin(alpha_0) + (y_obj - y_0).*cos(alpha_0));
% X_b = (x_obj - x_0)*H./(f*sin(alpha_0) + (y_obj - y_0).*cos(alpha_0));
% 
% Y_b = - Y_b;

% get closest point in the real framework 
idx_x = floor(x_obj*2/W);
idx_y = floor(y_obj*2/W);

signal = squeeze(dz(idx_x,idx_y,:));
end

function [TF,TF_pos,f] = fft_1D(signal,fs,padding_bool,add_pow2)
    % Computes the Fourier transform of a 1D vector
    % Inputs : 
    % - signal : 1D vector 
    % - fs : sampling frequency
    % - padding_bool : boolean to choose to pad or not
    % - add_pow2 : additional power of 2 if padding 
    
    % Outputs : 
    % - TF : Fourier transform 
    % - TF_pos : Fourier transform keeping only positive frequencies
    % - f : frequency array, only positive frequencies 
    
    original_length = numel(signal);
    padding_length = 2^(nextpow2(original_length) + add_pow2);

    signal = signal - mean(signal);
    if padding_bool 
        TF = fft(signal,padding_length);
        N = padding_length ;
        disp('Padding used')
    else 
        TF = fft(signal);
        N = original_length;
        disp('No padding')
    end

    TF = TF/original_length; % normalization of the FFT
    TF_pos = TF(1:N/2+1);
    TF_pos(2:end-1) = 2*TF_pos(2:end-1); % multiply by 2 for peaks that are both in positive an negative frequencies
    f = fs*(0:(N/2))/N; % frequency array
  
end 
