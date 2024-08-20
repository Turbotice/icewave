%%

%% Load matrix of signal from drones and buoys 
base = 'W:/SagWin2024/Data/0211/Drones/Noise_drone/';
filename = [base, 'Signal_initial_buoys_position.mat'];

disp('Loading data...')
load(filename)
disp('Data loaded')

%% Folder where to save images 
fig_folder = base; 

%% Filter each signal using a low pass filter 
fcut = 1.0; % cutoff frequency 
fs = B_matrix(2,2).fps; %acquisition frequency of buoys
filter_order = 6;
[b,a] = butter(filter_order,fcut/(fs/2));
filtered = filtfilt(b,a,B_matrix(3,2).signal);

[b,a] = butter(filter_order,0.6/(fs/2));
highfiltered = filtfilt(b,a,B_matrix(3,2).signal);

figure, 
% plot(B_matrix(3,2).UTC_t,B_matrix(3,2).signal,'LineWidth',2)
plot(B_matrix(3,2).UTC_t,highfiltered,'LineWidth',2)
hold on 
plot(B_matrix(3,2).UTC_t,filtered,'LineWidth',2)
grid on 
ax = gca;
ax.FontSize = 13;
set_Papermode(gcf)

%% Filter all signals 

fcut = 1.0; 
fs = B_matrix(2,2).fps; %acquisition frequency of buoys
filter_order = 6;
[b,a] = butter(filter_order,fcut/(fs/2));



for i = 1:3 % loop over buoys
    for j = 2:3 % loop over drones 
        B_matrix(i,j).filtered_signal = filtfilt(b,a,B_matrix(i,j).signal);
        B_matrix(i,j).filter.fcut = fcut;
        B_matrix(i,j).filter.coef = [b,a];
    end 
end 

%%
buoy = ['1','3','5'];

for i = 1:3

    figure(1),
    plot(B_matrix(i,2).UTC_t,B_matrix(i,2).filtered_signal,'LineWidth',1)
    hold on 
    plot(B_matrix(i,3).UTC_t,B_matrix(i,3).filtered_signal,'LineWidth',1)
    hold on 
    plot(B_matrix(i,1).UTC_t,B_matrix(i,1).filtered_signal,'LineWidth',1)
    grid on 
    xlabel('$\rm Time \: (UTC)$','Interpreter','latex')
    ylabel('$V_z \: \rm (m.s{-1})$','Interpreter','latex')
    legend('Bernache','Mesange')
    ylim([-1.0 1.0])
    ax = gca;
    ax.FontSize = 13;
    set_Papermode(gcf)

    figname = [fig_folder, 'Filtered_signal_buoy_', buoy(i)];
    saveas(gcf,figname,'fig')
    saveas(gcf,figname,'pdf')
    hold off
end 

%% Load buoys UTC_time 


base_buoys = '/media/turbots/DATA/thiou/labshared2/SagWin2024/Data/0211/BoueeVague/B1/mat/';

buoy_filename = [base_buoys 'buoy1_sbg_20240211_2000.mat'];

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

t_start = datetime(2024,2,11,20,34,50,'TimeZone','UTC'); 
t_end = datetime(2024,2,11,20,37,26,'TimeZone','UTC');

mask = (UTC_t >= t_start) & (UTC_t < t_end);
ROI_t = UTC_t(mask);

%%
for i = 1:3
    for j = 1:3
        B_matrix(i,j).UTC_t_buoy = ROI_t;
    end 
end 

%% Superpose filtered signal for all buoys

buoy = ['1','3','5'];

for i = 1:3

    figure(10),
    plot(B_matrix(i,2).UTC_t,B_matrix(i,2).filtered_signal,'LineWidth',1)
    hold on 
    plot(B_matrix(i,3).UTC_t,B_matrix(i,3).filtered_signal,'LineWidth',1)
    hold on 
    plot(B_matrix(i,1).UTC_t,B_matrix(i,1).signal,'k','LineWidth',1)
    grid on 
    xlabel('$\rm Time \: (UTC)$','Interpreter','latex')
    ylabel('$V_z \: \rm (m.s{-1})$','Interpreter','latex')
    legend('Bernache','Mesange','Buoy')
    ylim([-1.0 1.0])
    ax = gca;
    ax.FontSize = 13;
    set_Papermode(gcf)

    figname = [fig_folder, 'Filtered_buoy_superposition_signal_buoy_', buoy(i)];
    saveas(gcf,figname,'fig')
    saveas(gcf,figname,'pdf')
    hold off
end 


%% FFT over several pixels 

% load data from bernache 
%%

%# Buoy #1
% obj_txt = 'buoy_1';
% x_bernache = 1284;
% y_bernache = 752; 
% 
% x_mesange = 2570;
% y_mesange = 821;

%# Buoy #4
% obj_txt = 'buoy_4';
% x_bernache = 2593;
% y_bernache = 749;
% 
% x_mesange = 1126;
% y_mesange = 891;

%# Buoy #5
obj_txt = 'buoy_5';
x_bernache = 3204;
y_bernache = 759;

x_mesange = 586;
y_mesange = 911;

W = 32; % window size for PIV
window = 10; % semi width of window 
buoy_idx(3,:) = [floor(x_bernache*2/W),floor(y_bernache*2/W)]; % initial indices in field displacement matrix 

local_field = s_bernache.s.dz(buoy_idx(1,1)-window:buoy_idx(1,1)+window,buoy_idx(1,2)-window:buoy_idx(1,2)+window,:); 
padding_bool = 1;
add_pow2 = 1;
[FFT_t,TF_spectrum,f] = temporal_FFT(local_field,padding_bool,add_pow2,B_matrix(1,2).fps);

figure, 
loglog(f,TF_spectrum)
hold on 
loglog(B_matrix(1,2).f,abs(B_matrix(1,2).TF_pos))
grid on 
axis([5e-2 30 1e-5 1])
xlabel('$\rm f \: (Hz)$','Interpreter','latex')
% ylabel('$\langle \hat{V_z} \rangle _{x,y}(f) \: \rm (m.s{-1})$','Interpreter','latex')
ylabel('$\hat{s}(f) \: \rm (m.s{-1})$','Interpreter','latex')
legend('Space average','Single point')
ax = gca;
ax.FontSize = 13;
set_Papermode(gcf)

figname = [fig_folder, 'Space_average_FFT_buoy_', buoy(3), 'window_10_comparison_single' ];
saveas(gcf,figname,'fig')
saveas(gcf,figname,'pdf')