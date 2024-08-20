%% Load structure of buoys signal and positions 
path = 'F:/Rimouski_2024/Data/2024/0211/Drones/';
filename = [path 'buoys_drone_superposition_data.mat'];

load(filename)
disp('Data loaded')

%% Load Vz matrix for each drone (full structure)
path = 'F:/Rimouski_2024/Data/2024/0211/Drones/bernache/matData/structured_data/';
filename = [path 'PIV_processed_i011500_N15500_Dt4_b1_W32_full_scaled.mat'];

m_bernache = load(filename);

path = 'F:/Rimouski_2024/Data/2024/0211/Drones/mesange/matData/2-stereo_001/structured_data/';
filename = [path 'PIV_processed_i011496_N15496_Dt4_b1_W32_full_scaled.mat'];

m_mesange = load(filename);

%% Create a large structure combining both structures from bernache and mesange 

M = struct('bernache',m_bernache.m,'mesange',m_mesange.m);

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
buoy_key = 'B1'; % buoy to study
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

%% Average data using surrounding boxes 
W = 32; % bow width used for PIV (in pixels)
buoy_key = 'B1';
drone_key = 'bernache';
window = 2; % semi width of window 
frame = 1;
xpix = S.(buoy_key).(drone_key).pixel_pos(1);
ypix = S.(buoy_key).(drone_key).pixel_pos(2);
selected_boxes = [floor(xpix*2/W) - window : floor(xpix*2/W) + window;...
                  floor(ypix*2/W) - window : floor(ypix*2/W) + window]; 
ROI_X = M.(drone_key).X(selected_boxes(1,:),selected_boxes(2,:));
ROI_Y = M.(drone_key).Y(selected_boxes(1,:),selected_boxes(2,:));

figure, 
pcolor(M.(drone_key).X,M.(drone_key).Y,M.(drone_key).Vz(:,:,frame))
shading interp 
hold on 
plot(ROI_X,ROI_Y,'ro')
              

ROI_Vz = M.(drone_key).Vz(selected_boxes(1,:),selected_boxes(2,:),:);
signal_averaged = squeeze(mean(mean(ROI_Vz,2),1));

figure, 
plot(S.(buoy_key).(drone_key).UTC_t,S.(buoy_key).(drone_key).signal)
hold on 
plot(S.(buoy_key).(drone_key).UTC_t,signal_averaged)




%% Save new version of structure S (if needed) 







