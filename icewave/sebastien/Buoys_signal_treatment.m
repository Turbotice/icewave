% ###############################################################
%% ################## Load Buoys signal #########################
% ###############################################################

base_buoys = 'W:/SagWin2024/Data/0211/BoueeVague/B2/mat/';

buoy_filename = [base_buoys 'buoy2_sbg_20240211_1900.mat'];

disp('Loading data..')
B2 = load(buoy_filename);
disp('Data loaded')

buoy_filename = ['W:/SagWin2024/Data/0211/BoueeVague/B1/mat/' , 'buoy1_sbg_20240211_2000.mat'];
B1 = load(buoy_filename);
%% Plot raw data from buoys
IMU = B2.IMU;

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
plot(UTC_t,IMU.IMU_DATA.ACCEL_Z)
grid on 

%% Plot bernache and mesange signal with UTC_time 
%# Buoy #1
obj_txt = 'buoy_1';
x_bernache = 1284;
y_bernache = 752; 

x_mesange = 2570;
y_mesange = 821;

%# Buoy #2
% obj_tkt = 'buoy_2';
% x_bernache = 1457;
% y_bernache = 738;
% 
% x_mesange = 2401;
% y_mesange = 846;

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
vz = cumsum(filtered_buoy - mean(filtered_buoy));
figure, 
plot(ROI_t,vz/fs)

%% Superposition of drones and buoy with UTC time 
delta_t0 = s_bernache.s.UTC_time(1) - ROI_t(1); % difference of t0 between bernache UTC time and buoys UTC time 
add_dt = milliseconds(500); % add 500ms 
UTC_drone = s_bernache.s.UTC_time - delta_t0 + add_dt;

figure, 
plot(UTC_drone',signal_bernache)
hold on 
plot(UTC_drone',signal_mesange)
hold on 
plot(ROI_t,vz/fs,'k','LineWidth',1.0)
grid on 
xlabel('$t$','Interpreter','latex')
ylabel('$V_z \: \rm (m.s^{-1})$','Interpreter','latex')
title('Buoy 2','Interpreter','latex')
legend('Bernache','Mesange','Buoy','Interpreter','latex')
ylim([-1 1])
set_Papermode(gcf)
ax = gca;
ax.FontSize = 13;