%% Create a matrix of buoys position 
Buoy_pix = zeros([3,2,2]);
% dimensions : #1 : idx buoy {B1, B4, B5}, #2 : idx drone {bernache, mesange}
% #3 : (x,y) in pixels 

% Buoy #1
Buoy_pix(1,1,:) = [1284 , 752];
Buoy_pix(1,2,:) = [2570 , 821];

% Buoy #4
Buoy_pix(2,1,:) = [2593 , 749];
Buoy_pix(2,2,:) = [1126 , 891];

% Buoy #5 
Buoy_pix(3,1,:) = [3204 , 759];
Buoy_pix(3,2,:) = [586 , 911];

%% Load buoys matrix (synchronizing is not perfect) 
path = 'F:/Rimouski_2024/Data/2024/0211/Drones/' ;
filename = [path '1D_signal_drone_buoys_initial_positions.mat'];

load(filename)

%% Load Drone data 
path = 'F:/Rimouski_2024/Data/2024/0211/Drones/bernache/matData/structured_data/';
filename = [path 'PIV_processed_i011500_N15500_Dt4_b1_W32_full_scaled.mat'];

m_bernache = load(filename);

path = 'F:/Rimouski_2024/Data/2024/0211/Drones/mesange/matData/2-stereo_001/structured_data/';
filename = [path 'PIV_processed_i011496_N15496_Dt4_b1_W32_full_scaled.mat'];

m_mesange = load(filename);

%% Load single buoy data and extract synchronized vertical velocity 
key = 'B5';
table_buoy = struct('B1',1,'B4',2,'B5',3);
base = ['F:/Rimouski_2024/Data/0211/BoueeVague/' key '/mat/'];

buoy_filename = [base 'buoy' key(2) '_sbg_20240211_2000.mat'];

disp('Loading data..')
load(buoy_filename)
disp('Data loaded')

%%
figure(1), 
plot(IMU.IMU_DATA.ACCEL_Z)

%% Select the UTC-time associated to the beginning of the acquisition

Y = IMU.UTC_TIME.YEAR(1);
M = IMU.UTC_TIME.MONTH(1);
D = IMU.UTC_TIME.DAY(1);
H = IMU.UTC_TIME.HOUR(1);
MIN = IMU.UTC_TIME.MIN(1);
Sec = IMU.UTC_TIME.SEC(1);
% MS = IMU.UTC_TIME.NANOSEC(1);

% initial time of recording for buoys
t0 = datetime(Y,M,D,H,MIN,Sec,'TimeZone','UTC'); 
t0.TimeZone = 'UTC'; % converts time to UTC time 
t0.Format = 'yyyy-MMM-dd HH:mm:ss';

% Convert this UTC time to time since epoch 

t0_epoch = convertTo(t0,'posixtime');
% create an array of time t since the epoch for the acquisition
t_epoch = t0_epoch + IMU.IMU_DATA.t ; 

% create an array of UTC time for the acquisition 
UTC_t = datetime(t_epoch, 'convertfrom', 'posixtime', 'Format', 'MM/dd/yy HH:mm:ss.SSS');
UTC_t.TimeZone = 'UTC';

%% select relevant signal 
t_start = datetime(2024,2,11,20,34,50,'TimeZone','UTC'); 
t_end = datetime(2024,2,11,20,37,26,'TimeZone','UTC');

mask = (UTC_t >= t_start) & (UTC_t < t_end);
ROI_signal = IMU.IMU_DATA.ACCEL_Z(mask);
ROI_t = UTC_t(mask);

figure(2), 
plot(ROI_t,ROI_signal)

%% Applies a filter and Computes vertical velocity of boys
fs = 50; %acquisition frequency of buoys
fc = 0.1; %cutoff frequency 

[b,a] = butter(6,fc/(fs/2),"high");
% [b,a] = butter(6,/25,"high");
filtered_buoy = filtfilt(b,a,double(ROI_signal));

vz = cumsum(filtered_buoy - mean(filtered_buoy))/fs;
figure(3), 
plot(ROI_t,vz)
grid on 
%% Synchronize UTC of drone and buoys using buoy #5
% result : add_dt = 670 ms

delta_t0 = m_bernache.m.UTC_t(1) - ROI_t(1); % difference of t0 between bernache UTC time and buoys UTC time 
dt = (400:10:800);
add_dt = milliseconds(dt); % add 500ms 
key = 'B5';
t_min = datetime(2024,2,11,20,35,46,'TimeZone','UTC'); 
t_max = datetime(2024,2,11,20,36,10,'TimeZone','UTC');
signal_b = get_1D_signal(S.(key).bernache.pixel_pos(1),S.(key).bernache.pixel_pos(2),m_bernache.m.Vz,32);
signal_m = get_1D_signal(S.(key).mesange.pixel_pos(1),S.(key).mesange.pixel_pos(2),m_mesange.m.Vz,32);

for i = 1 : length(dt)
    UTC_drone(i,:) = m_bernache.m.UTC_t - delta_t0 + add_dt(i);
    
    mask_drone = (UTC_drone(i,:) >= t_min) & (UTC_drone(i,:) < t_max); 
    mask_buoy = (ROI_t >= t_min) & (ROI_t < t_max); 
    short_sig_b = signal_b(mask_drone);
    short_sig_m = signal_m(mask_drone);
    short_sig_buoy = vz(mask_buoy);
    
    I1 = interp1(ROI_t(mask_buoy),short_sig_buoy,UTC_drone(i,mask_drone));
    d(i,1) = sum(abs(I1 - short_sig_b'),'omitnan');
    d(i,2) = sum(abs(I1 - short_sig_m'),'omitnan');
end 

[~,idx_min_m] = min(d(:,2));
[~,idx_min_b] = min(d(:,1));

disp(add_dt(idx_min_m))
disp(add_dt(idx_min_b))
%% Modify UTC time for buoy, using bernache UTC as reference 

dt = 670;
add_dt = milliseconds(dt); % add 500ms 

t_min = datetime(2024,2,11,20,35,46,'TimeZone','UTC'); 
t_max = datetime(2024,2,11,20,36,10,'TimeZone','UTC');
signal_b = get_1D_signal(S.(key).bernache.pixel_pos(1),S.(key).bernache.pixel_pos(2),m_bernache.m.Vz,32);
signal_m = get_1D_signal(S.(key).mesange.pixel_pos(1),S.(key).mesange.pixel_pos(2),m_mesange.m.Vz,32);

UTC_drone = m_bernache.m.UTC_t - delta_t0 + add_dt;
UTC_buoy = ROI_t + delta_t0 - add_dt;

figure(4), 
plot(m_bernache.m.UTC_t,signal_b)
hold on 
plot(m_bernache.m.UTC_t,signal_m)
hold on 
plot(UTC_buoy,vz)
grid on 


%%
idx = table_buoy.(key);
B_matrix(idx,1) = struct('signal',vz,'UTC_t',UTC_buoy,'fps',50);
B_matrix(idx,2) = struct('signal',signal_b,'UTC_t',m_bernache.m.UTC_t,'fps',29.97);
B_matrix(idx,3) = struct('signal',signal_m,'UTC_t',m_bernache.m.UTC_t,'fps',29.97);

%% Save Buoy matrix containing signals 
base_save = 'F:/Rimouski_2024/Data/2024/0211/Drones/';
savename = [base_save '1D_signal_drone_buoys_initial_positions.mat'];

save(savename,'B_matrix','-v7.3')
disp('Data saved')


%% Create a matrix of signals 

S.B1.buoy.signal = B_matrix(1,1).signal;
S.B4.buoy.signal = B_matrix(2,1).signal;
S.B5.buoy.signal = B_matrix(3,1).signal;
W = 32;

buoys_idx = {'B1','B4','B5'};
drone_idx = {'bernache','mesange'};

for i = 1:3
    key = buoys_idx{i};
    S.(key).('buoy').('fps') = 50; % fps of buoys
    S.(key).('buoy').('UTC_t') = B_matrix(i,1).UTC_t; % UTC time of buoys 
    for j = 1:2
        drone_key = drone_idx{j};
        if j == 1
            Vz = m_bernache.m.Vz;
        else
            Vz = m_mesange.m.Vz;
        end    
        S.(key).(drone_key).('pixel_pos') = squeeze(Buoy_pix(i,j,:));
        S.(key).(drone_key).('fps') = 29.97;
        S.(key).(drone_key).('UTC_t') = B_matrix(i,j+1).UTC_t;
        S.(key).(drone_key).('signal') = get_1D_signal(Buoy_pix(i,j,1), Buoy_pix(i,j,2),Vz,W);
    end 
end 

%% Save structure 
base_save = 'F:/Rimouski_2024/Data/2024/0211/Drones/';
savename = [base_save 'buoys_drone_superposition_data.mat'];

save(savename,'S','-v7.3')
disp('Data saved')

%% Plot signals for a given buoy
buoy_key = 'B1';

figure, 
plot(S.(buoy_key).buoy.UTC_t,S.(buoy_key).buoy.signal,'k')
hold on 
plot(S.(buoy_key).bernache.UTC_t,S.(buoy_key).bernache.signal)
hold on 
plot(S.(buoy_key).mesange.UTC_t,S.(buoy_key).mesange.signal)

%% Apply a bandpass filter 

fcut = [0.1 , 1.0];
buoy_key = 'B1';
drone_key = 'bernache';
filtered = bandpass(S.(buoy_key).(drone_key).signal,fcut,S.(buoy_key).(drone_key).fps);

figure, 
plot(S.(buoy_key).(drone_key).UTC_t,filtered)
hold on 
plot(S.(buoy_key).(drone_key).UTC_t,S.(buoy_key).(drone_key).signal)
grid on 

%% Superpose filtered data and buoys signal 
filtered = struct('bernache',[],'mesange',[]);

for i  = 1:2
    drone_key = drone_idx{i};
    filtered.(drone_key) = bandpass(S.(buoy_key).(drone_key).signal,fcut,...
        S.(buoy_key).(drone_key).fps);
end 

figure, 
plot(S.(buoy_key).('bernache').UTC_t,filtered.('bernache'))
hold on 
plot(S.(buoy_key).('mesange').UTC_t,filtered.('mesange'))
hold on 
plot(S.(buoy_key).buoy.UTC_t,S.(buoy_key).buoy.signal,'k','LineWidth',1)
grid on 

%% Select time segment of Vz values for each drone 

s_bernache = m_bernache.m;
frame_max = 1800;

s_bernache.Vx = s_bernache.Vx(:,:,1:frame_max);
s_bernache.Vy = s_bernache.Vy(:,:,1:frame_max);
s_bernache.Vz = s_bernache.Vz(:,:,1:frame_max);
s_bernache.t = s_bernache.t(1:frame_max);
s_bernache.UTC_t = s_bernache.UTC_t(1:frame_max);

base_save = 'F:/Rimouski_2024/Data/2024/0211/Drones/bernache/matData/structured_data/';
filename = [base_save 'PIV_data_scaled_bernache_small_size.mat'];
save(filename, 's_bernache','-v7.3')
disp('Reduced structure saved')

%%
s_mesange = m_mesange.m;
frame_max = 1800;

s_mesange.Vx = s_mesange.Vx(:,:,1:frame_max);
s_mesange.Vy = s_mesange.Vy(:,:,1:frame_max);
s_mesange.Vz = s_mesange.Vz(:,:,1:frame_max);
s_mesange.t = s_mesange.t(1:frame_max);
s_mesange.UTC_t = s_mesange.UTC_t(1:frame_max);

base_save = 'F:/Rimouski_2024/Data/2024/0211/Drones/mesange/matData/2-stereo_001/structured_data/';
filename = [base_save 'PIV_data_scaled_mesange_small_size.mat'];
save(filename, 's_mesange','-v7.3')
disp('Reduced structure saved')

% figure,
% pcolor(m.X,m.Y,m.Vz(:,:,idx_frame))
% shading interp






% #############################
%% FUNCTION SECTION 
% #############################

function [signal] = get_1D_signal(x_obj,y_obj,Vz,W)

% This function computes the vertical velocity at a given position of the
% PIV field, and returns it as a 1D signal
% The function takes the following arguments : 
% - x_obj : x-coordinate on the camera (in pixels)
% - y_obj : y-coordinate on the camera (in pixels)
% - Vz : vertical velocity matrix [nx,ny,nt]
% - W : window size used for the PIV processing (pixels)

% object position in the sea ice framework
% Y_b = (y_obj - y_0)*H/sin(alpha_0)./(f*sin(alpha_0) + (y_obj - y_0).*cos(alpha_0));
% X_b = (x_obj - x_0)*H./(f*sin(alpha_0) + (y_obj - y_0).*cos(alpha_0));
% 
% Y_b = - Y_b;

% get closest point in the real framework 
idx_x = floor(x_obj*2/W);
idx_y = floor(y_obj*2/W);

signal = squeeze(Vz(idx_x,idx_y,:));
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




