clear all 
close all
%% This script describes a method to extract objects position of a given
% colour in a 2D field taken by drone

% ########################
% EXTRACTION FOR FULMAR 
% ########################

% Load the whole video 
base = 'W:/SagWin2024/Data/0211/Drones/Fulmar/';
filename = fullfile(base , 'SWO_FUL_20240211T203233UTC.mp4');
% filename = 'LEGER_SWO_FUL_20240211T203233UTC.MP4'
v = VideoReader(filename);

%% Save data in a .mat file 

% selected values on each channel 
R0 = 0.83; % value on canal R
G0 = 0.34; % value on canal G 
B0 = 0.31; % value on canal B 

RGB0 = [R0,G0,B0];
% threshold for binarization 
threshold = 0.8;

% Scaling factor 
h_drone = 140; % drone altitude in meter
theta_x = 32.75; % AFOV of drone in (°)
L_x = 3840; % number of pixels along x-direction

fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
fps = v.FrameRate;
% define a structuring element for erosion 
SE = strel('square',2);
min_area = 4;
min_distance = 20; % minimal distance in pixels to merge objects
idx_first_frame = 4500; % index of first frame used for PIV
Nb_frames = v.NumFrames - idx_first_frame + 1;
Nan_array = NaN(Nb_frames,1); % create an array full of Nan values


base_save = [base 'buoys_tracking_0211_V2/'];
if ~exist(base_save)
    mkdir(base_save)
end 

file2save = [base_save 'Buoys_tracking_pix_fulmar.mat'];
heads = {'idx','t','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','x6','y6'};
% create an empty structure based on heads
S_buoys = struct();

for i = 1:length(heads)
    S_buoys.(heads{i}) = Nan_array;
end

selected_frames = (idx_first_frame : 1 : v.NumFrames);

%% Process each image

for i = 1:Nb_frames

    i0 = selected_frames(i);
    idx_frame = i0; 
    disp(i0)
    % Extract frame i from the video 
    img = read(v,i0);
%     img = readFrame(v);

    new_table = detect_buoys(img,RGB0,threshold,SE,min_distance);
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    
    t = (i0 - 1)/fps; % time 

    S_buoys.idx(i) = idx_frame;
    S_buoys.t(i) = t;


    for j = 1:length(new_table.Area)
        if new_table.Area(j) > min_area
            if j == 1
                count_obj = 1;
            else
                count_obj = count_obj + 1;
            end 
            xc = new_table.Centroid(j,1);
            yc = new_table.Centroid(j,2);
            head_xc = heads{2*count_obj + 1};
            S_buoys.(head_xc)(i) = xc;
            head_yc = heads{2*count_obj + 2};
            S_buoys.(head_yc)(i) = yc;
        end
    end 

end 
save(file2save,'S_buoys','-v7.3')
disp('Done.')

% #############################
%% Extraction for bernache 
% #############################

% Load the whole video 
base = 'F:/Rimouski_2024/Data/0211/Drones/bernache/stereo_001/';
addpath('Z:/skuchly/Codes/Drone/Function_package/');

filename = fullfile(base , 'DJI_20240211153234_0273_D.MP4');
% filename = 'LEGER_SWO_FUL_20240211T203233UTC.MP4'
v = VideoReader(filename);

filename = fullfile(base , 'DJI_20240211153621_0274_D.MP4');
v2 = VideoReader(filename);

% Load flight_parameters 
base_param = 'W:/SagWin2024/Data/0211/Drones/bernache/matData/18-stereo_001/' ;
param_file = [base_param 'Param_Debut_bernache_PIV_DJI_20240211153234_0273_D.mat'];
load(param_file)
%%
% define virtual frame 
first_frame = 11500; 
last_frame = 6803*2 + v2.NumFrames;
virtual_frames = (first_frame : last_frame);

change_video = 6803*2; % last frame of video #2
total_frame = round(v.Duration*v.FrameRate);

v0_NumFrames = 6803;
% selected values on each channel 
R0 = 0.83; % value on canal R
G0 = 0.34; % value on canal G 
B0 = 0.31; % value on canal B 

RGB0 = [R0,G0,B0];

alpha_0 = param.alpha_0*pi/180;
f = 2700; % focale in pixel 
img_0 = read(v,1);
x_0 = (size(img_0,2) + 1)/2;
y_0 = (size(img_0,1) + 1)/2;

% threshold for binarization 
threshold = 0.8;

fps = v.FrameRate;
% define a structuring element for erosion 
SE = strel('square',2);
min_area = 3; % minimal area of detected objects 
min_distance = 20; % minimal distance in pixels to merge objects

base_save = 'W:/SagWin2024/Data/0211/Drones/bernache/buoys_tracking_0211_V2/';
if ~exist(base_save)
    mkdir(base_save)
end 

Nb_frames = length(virtual_frames);
Nan_array = NaN(Nb_frames,1); % create an array full of Nan values

file2save = [base_save 'Buoys_tracking_pix_bernache.mat'];
heads = {'idx','t','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','x6','y6'};
% create an empty structure based on heads
S_buoys = struct();

for i = 1:length(heads)
    S_buoys.(heads{i}) = Nan_array;
end
%% Loop over all frames 
for i0 = 1:length(virtual_frames)
    virtual_index = virtual_frames(i0); 
    
    if virtual_index <= change_video
        i = virtual_index - v0_NumFrames;
        % Extract frame i from the video 
        img = read(v,i);
        disp(i)
        
    else 
        i = virtual_index - change_video;
        img = read(v2,i);
        disp(i)
    end 

    new_table = detect_buoys(img,RGB0,threshold,SE,min_distance);
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    t = (i0 - 1)/fps; % time 

    S_buoys.idx(i0) = virtual_index;
    S_buoys.t(i0) = t;
    
    for j = 1:length(new_table.Area)
        if new_table.Area(j) > min_area
            if j == 1
                count_obj = 1;
            else
                count_obj = count_obj + 1;
            end 
            xc = new_table.Centroid(j,1);
            yc = new_table.Centroid(j,2);
            head_xc = heads{2*count_obj + 1};
            S_buoys.(head_xc)(i0) = xc;
            head_yc = heads{2*count_obj + 2};
            S_buoys.(head_yc)(i0) = yc;
        end
    end 
    
end 
save(file2save,'S_buoys','-v7.3')
disp('Done.')

% #############################
%% Extraction for mesange 
% #############################

% Load the whole video 
base = 'F:/Rimouski_2024/Data/0211/Drones/mesange/stereo_001/';
addpath('Z:/skuchly/Codes/Drone/Function_package/');

filename = fullfile(base , 'DJI_20240211213233_0222_D.MP4');
% filename = 'LEGER_SWO_FUL_20240211T203233UTC.MP4'
v = VideoReader(filename);

filename = fullfile(base , 'DJI_20240211213620_0223_D.MP4');
v2 = VideoReader(filename);

filename = fullfile(base , 'DJI_20240211212846_0221_D.MP4');
v0 = VideoReader(filename);

% Load flight_parameters 
base_param = 'W:/SagWin2024/Data/0211/Drones/mesange/matData/2-stereo_001/' ;
param_file = [base_param 'Param_Debut_mesange_PIV_DJI_20240211213233_0222_D.mat'];
load(param_file)

%%
% define virtual frame 
first_frame = 11496; 
last_frame = v0.NumFrames + v.NumFrames + v2.NumFrames;
virtual_frames = (first_frame : last_frame);

change_video = v0.NumFrames + v.NumFrames; % last frame of video #2
total_frame = round(v.Duration*v.FrameRate);

v0_NumFrames = v0.NumFrames;
% selected values on each channel 
R0 = 0.83; % value on canal R
G0 = 0.34; % value on canal G 
B0 = 0.31; % value on canal B 

alpha_0 = param.alpha_0*pi/180;
f = 2700; % focale in pixel 
img_0 = read(v,1);
x_0 = (size(img_0,2) + 1)/2;
y_0 = (size(img_0,1) + 1)/2;

% threshold for binarization 
threshold = 0.8;

fps = v.FrameRate;
% define a structuring element for erosion 
SE = strel('square',2);
min_area = 3;
min_distance = 20; % minimal distance in pixels to merge objects

base_save = 'W:/SagWin2024/Data/0211/Drones/mesange/buoys_tracking_0211_V2/';
if ~exist(base_save)
    mkdir(base_save)
end 

Nb_frames = length(virtual_frames);
Nan_array = NaN(Nb_frames,1); % create an array full of Nan values

file2save = [base_save 'Buoys_tracking_pix_mesange.mat'];
heads = {'idx','t','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','x6','y6'};
% create an empty structure based on heads
S_buoys = struct();

for i = 1:length(heads)
    S_buoys.(heads{i}) = Nan_array;
end

%%
for i0 = 1:length(virtual_frames)
    virtual_index = virtual_frames(i0); 
    
    if virtual_index <= change_video
        i = virtual_index - v0_NumFrames;
        % Extract frame i from the video 
        img = read(v,i);
        disp(i)
        
    else 
        i = virtual_index - change_video;
        img = read(v2,i);
        disp(i)
    end 

    new_table = detect_buoys(img,RGB0,threshold,SE,min_distance);
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    
    t = (i0 - 1)/fps; % time 
    S_buoys.idx(i0) = virtual_index;
    S_buoys.t(i0) = t;
    
    for j = 1:length(new_table.Area)
        if new_table.Area(j) > min_area
            if j == 1
                count_obj = 1;
            else
                count_obj = count_obj + 1;
            end 
            xc = new_table.Centroid(j,1);
            yc = new_table.Centroid(j,2);
            head_xc = heads{2*count_obj + 1};
            S_buoys.(head_xc)(i0) = xc;
            head_yc = heads{2*count_obj + 2};
            S_buoys.(head_yc)(i0) = yc;
        end
    end 
    
end 
save(file2save,'S_buoys','-v7.3')
disp('Done.')











%#########################################################################################################
%#################################### OLD VERSION ######################################################################
%##########################################################################################################
%##########################################################################################################
%##########################################################################################################

%%
% selected values on each channel 
R0 = 0.83; % value on canal R
G0 = 0.34; % value on canal G 
B0 = 0.31; % value on canal B 

RGB0 = [R0,G0,B0];
% threshold for binarization 
threshold = 0.8;

% Scaling factor 
h_drone = 140; % drone altitude in meter
theta_x = 32.75; % AFOV of drone in (°)
L_x = 3840; % number of pixels along x-direction

fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
fps = v.FrameRate;
% define a structuring element for erosion 
SE = strel('square',2);
min_area = 1;
min_distance = 20; % minimal distance in pixels to merge objects

base_save = [base 'buoys_tracking_0211_V2/'];
if ~exist(base_save)
    mkdir(base_save)
end 
% create a csv file and write heads 
csv_file = [base_save 'Buoys_tracking_fulmar.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,','X6,','Y6,'];
writematrix(heads,csv_file,'WriteMode','append')

% create a csv file for positions in pixels 
csv_file_pix = [base_save 'Buoys_tracking_pix_fulmar.csv'];
heads = ['idx,','t,','x1,','y1,','x2,','y2,','x3,','y3,','x4,','y4,','x5,','y5,','x6,','y6,'];
writematrix(heads,csv_file_pix,'WriteMode','append')

%%

for i = 4500:v.NumFrames

    disp(i)
    % Extract frame i from the video 
    img = read(v,i);
%     img = readFrame(v);

    new_table = detect_buoys(img,RGB0,threshold,SE,min_distance);
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    
    t = (i - 1)/fps; % time 
    row = [i, t]; % initialize row to write in a csv file 
    row_pix = [i,t];
    for j = 1:length(new_table.Area)
        if new_table.Area(j) > min_area 
            xc = new_table.Centroid(j,1);
            yc = new_table.Centroid(j,2);
            % convert coordinates in meters 
            X_buoy = xc/fx_pix;
            Y_buoy = yc/fx_pix; 
            row = cat(2,row,X_buoy,Y_buoy);
            row_pix = cat(2,row_pix,xc,yc);
        end
    end 
    
    row = string(row);
    row_pix = string(row_pix);
    % Write positions of detected ROIs in a csv file 
    writematrix(row,csv_file,'WriteMode','append')
    writematrix(row_pix,csv_file_pix,'WriteMode','append')

end 
disp('Done.')

%% Save parameters 

file_param = [base_save 'Parameters_buoy_tracking_fulmar_0211'];
save(file_param,'RGB0','threshold','SE','min_area','min_distance','-v7.3')

%% Correction to the image center for Fulmar 
% Real frame work 

% load csv file 
path = 'W:/SagWin2024/Data/0211/Drones/Fulmar/buoys_tracking_V2/';
initial_file = [path 'Buoys_tracking_fulmar.csv'];

table = readtable(initial_file);

%%
% Scaling factor 
h_drone = 140; % drone altitude in meter
theta_x = 32.75; % AFOV of drone in (°)
L_x = 3840; % number of pixels along x-direction

fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter

X = table{:,3:14};

X0 = (3840 + 1)/2/fx_pix;
Y0 = (2160 + 1)/2/fx_pix;

% put center of the coordinate system on X0,Y0
X(:,1:2:11) = X(:,1:2:11) - X0;
X(:,2:2:12) = -(X(:,2:2:12) - Y0); % Y is oriented upward

% merge buoys positions with frame and time columns
Frame_time = table{:,1:2};

new_array = cat(2,Frame_time,X);
disp('New array computed')
%% Write new array in a new csv file 
% create a csv file and write heads 
csv_file = [path 'Buoys_tracking_center.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,','X6,','Y6,'];
writematrix(heads,csv_file,'WriteMode','append')

for i = 1 : size(new_array,1)
    disp(i)
    row = new_array(i,:);
    row = string(row);
%     Write positions of detected ROIs in a csv file 
    writematrix(row,csv_file,'WriteMode','append')
end
% #############################
%% Extraction for bernache 
% #############################

% Load the whole video 
base = 'F:/Rimouski_2024/Data/0211/Drones/bernache/stereo_001/';
addpath('Z:/skuchly/Codes/Drone/Function_package/');

filename = fullfile(base , 'DJI_20240211153234_0273_D.MP4');
% filename = 'LEGER_SWO_FUL_20240211T203233UTC.MP4'
v = VideoReader(filename);

filename = fullfile(base , 'DJI_20240211153621_0274_D.MP4');
v2 = VideoReader(filename);

% Load flight_parameters 
base_param = 'W:/SagWin2024/Data/0211/Drones/bernache/matData/18-stereo_001/' ;
param_file = [base_param 'Param_Debut_bernache_PIV_DJI_20240211153234_0273_D.mat'];
load(param_file)
%%
% define virtual frame 
first_frame = 11500; 
last_frame = 6803*2 + v2.NumFrames;
virtual_frames = (first_frame : last_frame);

change_video = 6803*2; % last frame of video #2
total_frame = round(v.Duration*v.FrameRate);

v0_NumFrames = 6803;
% selected values on each channel 
R0 = 0.83; % value on canal R
G0 = 0.34; % value on canal G 
B0 = 0.31; % value on canal B 

RGB0 = [R0,G0,B0];

alpha_0 = param.alpha_0*pi/180;
f = 2700; % focale in pixel 
img_0 = read(v,1);
x_0 = (size(img_0,2) + 1)/2;
y_0 = (size(img_0,1) + 1)/2;

% threshold for binarization 
threshold = 0.8;

fps = v.FrameRate;
% define a structuring element for erosion 
SE = strel('square',2);
min_area = 3; % minimal area of detected objects 
min_distance = 20; % minimal distance in pixels to merge objects

base_save = 'W:/SagWin2024/Data/0211/Drones/bernache/buoys_tracking_0211_V2/';
if ~exist(base_save)
    mkdir(base_save)
end 

% create a csv file and write heads 
csv_file = [base_save 'Buoys_tracking_bernache.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,','X6,','Y6,'];
writematrix(heads,csv_file,'WriteMode','append')

% create a csv file for positions in pixels 
csv_file_pix = [base_save 'Buoys_tracking_pix_bernache.csv'];
heads = ['idx,','t,','x1,','y1,','x2,','y2,','x3,','y3,','x4,','y4,','x5,','y5,','x6,','y6,'];
writematrix(heads,csv_file_pix,'WriteMode','append')

%% Loop over all frames 
for i0 = 1:length(virtual_frames)
    virtual_index = virtual_frames(i0); 
    
    if virtual_index <= change_video
        i = virtual_index - v0_NumFrames;
        % Extract frame i from the video 
        img = read(v,i);
        disp(i)
        
    else 
        i = virtual_index - change_video;
        img = read(v2,i);
        disp(i)
    end 

    new_table = detect_buoys(img,RGB0,threshold,SE,min_distance);
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    
    t = (i0 - 1)/fps; % time 
    row = [i, t]; % initialize row to write in a csv file 
    row_pix = [i,t];
    for j = 1:length(new_table.Area)
        if new_table.Area(j) > min_area 
            xc = new_table.Centroid(j,1);
            yc = new_table.Centroid(j,2);
            % convert coordinates in meters 
            
            [X_buoy,Y_buoy] = projection_real_space(xc,yc,x_0,y_0,param.H,alpha_0,f);
            
            row = cat(2,row,X_buoy,Y_buoy);
            row_pix = cat(2,row_pix,xc,yc);
        end
    end 
    
    row = string(row);
    row_pix = string(row_pix);
%     Write positions of detected ROIs in a csv file 
    writematrix(row,csv_file,'WriteMode','append')
    writematrix(row_pix,csv_file_pix,'WriteMode','append')

end 
disp('Done.')


%% Save parameters
file_param = [base 'Parameters_buoy_tracking_bernache_0211'];
save(file_param,'RGB0','threshold','SE','min_area','min_distance','-v7.3')

% #############################
%% Extraction for mesange 
% #############################

% Load the whole video 
base = 'F:/Rimouski_2024/Data/0211/Drones/mesange/stereo_001/';
addpath('Z:/skuchly/Codes/Drone/Function_package/');

filename = fullfile(base , 'DJI_20240211213233_0222_D.MP4');
% filename = 'LEGER_SWO_FUL_20240211T203233UTC.MP4'
v = VideoReader(filename);

filename = fullfile(base , 'DJI_20240211213620_0223_D.MP4');
v2 = VideoReader(filename);

filename = fullfile(base , 'DJI_20240211212846_0221_D.MP4');
v0 = VideoReader(filename);

% Load flight_parameters 
base_param = 'W:/SagWin2024/Data/0211/Drones/mesange/matData/2-stereo_001/' ;
param_file = [base_param 'Param_Debut_mesange_PIV_DJI_20240211213233_0222_D.mat'];
load(param_file)

%%
% define virtual frame 
first_frame = 11496; 
last_frame = v0.NumFrames + v.NumFrames + v2.NumFrames;
virtual_frames = (first_frame : last_frame);

change_video = v0.NumFrames + v.NumFrames; % last frame of video #2
total_frame = round(v.Duration*v.FrameRate);

v0_NumFrames = v0.NumFrames;
% selected values on each channel 
R0 = 0.83; % value on canal R
G0 = 0.34; % value on canal G 
B0 = 0.31; % value on canal B 

alpha_0 = param.alpha_0*pi/180;
f = 2700; % focale in pixel 
img_0 = read(v,1);
x_0 = (size(img_0,2) + 1)/2;
y_0 = (size(img_0,1) + 1)/2;

% threshold for binarization 
threshold = 0.8;

fps = v.FrameRate;
% define a structuring element for erosion 
SE = strel('square',2);
min_area = 3;
min_distance = 20; % minimal distance in pixels to merge objects

base_save = 'W:/SagWin2024/Data/0211/Drones/mesange/buoys_tracking_0211_V2/';
if ~exist(base_save)
    mkdir(base_save)
end 

% create a csv file and write heads 
csv_file = [base_save 'Buoys_tracking_mesange.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,'];
writematrix(heads,csv_file,'WriteMode','append')

% create a csv file for positions in pixels 
csv_file_pix = [base_save 'Buoys_tracking_pix_mesange.csv'];
heads = ['idx,','t,','x1,','y1,','x2,','y2,','x3,','y3,','x4,','y4,','x5,','y5,'];
writematrix(heads,csv_file_pix,'WriteMode','append')

%%
for i0 = 1:length(virtual_frames)
    virtual_index = virtual_frames(i0); 
    
    if virtual_index <= change_video
        i = virtual_index - v0_NumFrames;
        % Extract frame i from the video 
        img = read(v,i);
        disp(i)
        
    else 
        i = virtual_index - change_video;
        img = read(v2,i);
        disp(i)
    end 

    new_table = detect_buoys(img,RGB0,threshold,SE,min_distance);
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    
    t = (i - 1)/fps; % time 
    row = [i, t]; % initialize row to write in a csv file 
    row_pix = [i,t];
    for j = 1:length(new_table.Area)
        if (new_table.Area(j) > min_area) && (new_table.Centroid(j,1) > 200)
            xc = new_table.Centroid(j,1);
            yc = new_table.Centroid(j,2);
            % convert coordinates in meters 
            
            [X_buoy,Y_buoy] = projection_real_space(xc,yc,x_0,y_0,param.H,alpha_0,f);
            
            row = cat(2,row,X_buoy,Y_buoy);
            row_pix = cat(2,row_pix,xc,yc);
        end
    end 
    
    row = string(row);
    row_pix = string(row_pix);
%     Write positions of detected ROIs in a csv file 
    writematrix(row,csv_file,'WriteMode','append')
    writematrix(row_pix,csv_file_pix,'WriteMode','append')

end 
disp('Done.')

%% Save parameters
file_param = [base_save 'Parameters_buoy_tracking_mesange_0211'];
save(file_param,'RGB0','threshold','SE','min_area','min_distance','-v7.3')

% ##############################################
%% Modification of frame index and time arrays
% ##############################################

base = 'W:/SagWin2024/Data/0211/Drones/';
paths = {'Fulmar/buoys_tracking_0211_V2/','mesange/buoys_tracking_0211_V2/','bernache/buoys_tracking_0211_V2/'};
names = {'Buoys_tracking_pix_fulmar','Buoys_tracking_pix_mesange','Buoys_tracking_pix_bernache'};
indices_first_frame = [4500, 11496, 11500];

for i0 = 1:3
    filename = [base paths{i0} names{i0} '.csv'];
    % load table 
    table = readtable(filename);
    if i0 == 1 
        X = table{:,3:end-2};
    else
        X = table{:,3:end};
    end


    idx_frame = (indices_first_frame(i0) : indices_first_frame(i0) + size(X,1) -1)';

    t = (0:size(X,1)-1)';
    new_array = cat(2,idx_frame,t);
    new_array = cat(2,new_array,X);

    % create a csv file and write heads 
    csv_file = [base paths{i0} names{i0} '_pix_new.csv'];
    heads = ['idx_frame,','t,','x1,','y1,','x2,','y2,','x3,','y3,','x4,','y4,','x5,','y5,','x6,','y6,'];
    writematrix(heads,csv_file,'WriteMode','append')

    for i = 1 : size(X,1)
        disp(i)
        row = new_array(i,:);
        row = string(row);
    %     Write positions of detected ROIs in a csv file 
        writematrix(row,csv_file,'WriteMode','append')
    end
end 