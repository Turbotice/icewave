
clear all 
close all
%% This script describes a method to extract objects position of a given
% colour in a 2D field taken by drone

% Load the whole video 
base = 'W:/SagWin2024/Data/0211/Drones/Fulmar/';
% 'E:/Rimouski_2024/Data/2024/0211/Drones/Fulmar/FULMAR_vertical/';
% addpath(base);

filename = fullfile(base , 'SWO_FUL_20240211T203233UTC.mp4');
% filename = 'LEGER_SWO_FUL_20240211T203233UTC.MP4'
v = VideoReader(filename);
%%
i = 1000;
total_frame = round(v.Duration*v.FrameRate);

% selected values on each channel 
R0 = 0.83; % value on canal R
G0 = 0.34; % value on canal G 
B0 = 0.31; % value on canal B 

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
% create a csv file and write heads 
csv_file = [base 'Buoys_tracking.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,','X6,','Y6,'];
writematrix(heads,csv_file,'WriteMode','append')

% create a csv file for positions in pixels 
csv_file_pix = [base 'Buoys_tracking_pix.csv'];
heads = ['idx,','t,','x1,','y1,','x2,','y2,','x3,','y3,','x4,','y4,','x5,','y5,','x6,','y6,'];
writematrix(heads,csv_file_pix,'WriteMode','append')
%%

for i = 1:3

    disp(i)
    % Extract frame i from the video 
    img = read(v,i);
%     img = readFrame(v);

    % Build a distance from colours

    img_double = im2double(img);
    % build intensity matrix 
    I = sqrt((img_double(:,:,1) - R0).^2 + (img_double(:,:,2) - G0).^2 + (img_double(:,:,3) - B0).^2);
    % Mask large values 
    Inew = ones(size(I)) - I;
%     figure, 
%     imagesc(Inew)
%     colorbar()

    Inew(Inew >= threshold) = 1;
    Inew(Inew < threshold) = 0;

    % Detect ROI and extract there position 
    
%     CC = bwconncomp(Inew);
%     disp(['Number of detected objects :' num2str(CC.NumObjects)])
    Inew = imerode(Inew,SE); % erosion 
    CC = bwconncomp(Inew);

    stats = regionprops("table",CC);
% 
    figure, 
    imshow(img)
    hold on 
    plot(stats.Centroid(:,1),stats.Centroid(:,2),'or','MarkerSize',6)
    disp(stats)
 
%     stats.Centroid(:,1) = stats.Centroid(:,1);
%     stats.Centroid(:,2) = stats.Centroid(:,2);
    stats = sortrows(stats,'Centroid','descend'); % sort rows according to x position
    
    
    % Compute distance matrix between objects 
    D = squareform(pdist(stats.Centroid)); 
    
    % Merge particles that are too close to each other 
    mask_dist = (D < min_distance) & (D > 0);
    
    if ismember(1,mask_dist)
        [r,c] = find(mask_dist); % find indices of non-zero elements
        % select only superior part of the matrix 
        list_i = r(r>c);
        list_j = c(r>c);

        particle_i = list_i(~ ismember(list_i,list_j));
        particle_j = list_j(~ ismember(list_i,list_j));

        unique_part = unique(particle_i);

        new_part_centroid = zeros(length(unique_part),2);
        new_part_area = zeros(length(unique_part),1);
        for i0 = 1:length(unique_part)
            part = unique_part(i0);
            associated_part = particle_j(particle_i == part);
            list_particle = cat(1,part,associated_part);

            % computes mass center coordinates 
            new_part_centroid(i0,:) = sum(stats.Area(list_particle).*stats.Centroid(list_particle,:))...
                /sum(stats.Area(list_particle));
            new_part_area(i0) = sum(stats.Area(list_particle));
        end 
        % Get particles that do not need to be merged

        idx_list = (1:size(D,1));
        non_merged_mask = ~ (ismember(idx_list,particle_i) + ismember(idx_list,particle_j));
        non_merged = idx_list(non_merged_mask);

        non_merged_centroid = stats.Centroid(non_merged,:);
        non_merged_area = stats.Area(non_merged);
        new_centroid = cat(1,non_merged_centroid,new_part_centroid);
        new_area = cat(1,non_merged_area,new_part_area);
        
        Centroid = new_centroid;
        Area = new_area;
        new_table = table(Centroid,Area);

    else 
        
        new_table = stats;
    end 
    
    new_table = sortrows(new_table,'Centroid','descend');
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    
    %
%     for i0 = 1:length(list_i)
%         idx_i = list_i(i0);
%         idx_j = list_j(i0);
% 
%         new_table(i0,:) = (stats.Centroid(idx_i,:) + stats.Centroid(idx_j,:))/2;
%     end 
%     
%     mask = stats.Area > min_area;
%     Nb_objects = numel(stats.Area(mask));
%     disp(['Number of detected objects,after erosion :' num2str(Nb_objects)])
%     
    % get position scaled in meter 
%     X_buoys = zeros(Nb_objects,1);
%     Y_buoys = zeros(Nb_objects,1);
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
%     Write positions of detected ROIs in a csv file 
    writematrix(row,csv_file,'WriteMode','append')
    writematrix(row_pix,csv_file_pix,'WriteMode','append')

end 
disp('Done.')

%% Save parameters 

file_param = [base 'Parameters_buoy_tracking_fulmar_0211'];
save(file_param,'R0','B0','G0','threshold','SE','min_area','min_distance','-v7.3')

% #############################
%% Extraction for bernache 
% #############################

% Load the whole video 
base = 'E:/Rimouski_2024/Data/2024/0211/Drones/bernache/stereo_001/';
% 'E:/Rimouski_2024/Data/2024/0211/Drones/Fulmar/FULMAR_vertical/';
% addpath(base);

filename = fullfile(base , 'DJI_20240211153234_0273_D.MP4');
% filename = 'LEGER_SWO_FUL_20240211T203233UTC.MP4'
v = VideoReader(filename);

filename = fullfile(base , 'DJI_20240211153621_0274_D.MP4');
v2 = VideoReader(filename);

filename = fullfile(base , 'DJI_20240211152847_0272_D.MP4');
v0 = VideoReader(filename);

% Load flight_parameters 
base_param = 'W:/SagWin2024/Data/0211/Drones/bernache/matData/18-stereo_001/' ;
param_file = [base_param 'Param_Debut_bernache_PIV_DJI_20240211153234_0273_D.mat'];
load(param_file)
%%
% define virtual frame 
first_frame = 11500; 
last_frame = v0.NumFrames + v.NumFrames + v2.NumFrames;
virtual_frames = (first_frame : last_frame);

change_video = v0.NumFRames + v.NumFrames; % last frame of video #2
total_frame = round(v.Duration*v.FrameRate);

v0_NumFrames = v0.NumFrames;
% selected values on each channel 
R0 = 0.83; % value on canal R
G0 = 0.34; % value on canal G 
B0 = 0.31; % value on canal B 

alpha_0 = param.alpha_0*pi/180;
f = 2830; % focale in pixel 
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
% create a csv file and write heads 
csv_file = [base 'Buoys_tracking.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,','X6,','Y6,'];
writematrix(heads,csv_file,'WriteMode','append')

% create a csv file for positions in pixels 
csv_file_pix = [base 'Buoys_tracking_pix.csv'];
heads = ['idx,','t,','x1,','y1,','x2,','y2,','x3,','y3,','x4,','y4,','x5,','y5,','x6,','y6,'];
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

%     img = readFrame(v);

    % Build a distance from colours

    img_double = im2double(img);
    % build intensity matrix 
    I = sqrt((img_double(:,:,1) - R0).^2 + (img_double(:,:,2) - G0).^2 + (img_double(:,:,3) - B0).^2);

    % Mask large values 
    Inew = ones(size(I)) - I;

    Inew(Inew >= threshold) = 1;
    Inew(Inew < threshold) = 0;

    % Detect ROI and extract there position 
    
%     CC = bwconncomp(Inew);
%     disp(['Number of detected objects :' num2str(CC.NumObjects)])
    Inew = imerode(Inew,SE); % erosion 
    CC = bwconncomp(Inew);

    stats = regionprops("table",CC);
% 
%     figure, 
%     imagesc(Inew)
%     hold on 
%     plot(stats.Centroid(:,1),stats.Centroid(:,2),'or')
%     disp(stats)
 
%     stats.Centroid(:,1) = stats.Centroid(:,1);
%     stats.Centroid(:,2) = stats.Centroid(:,2);
    stats = sortrows(stats,'Centroid','ascend'); % sort rows according to x position
    
    
    % Compute distance matrix between objects 
    D = squareform(pdist(stats.Centroid)); 
   
    % Merge particles that are too close to each other 
    mask_dist = (D < min_distance) & (D > 0);
    
    if ismember(1,mask_dist)
        [r,c] = find(mask_dist); % find indices of non-zero elements
        % select only superior part of the matrix 
        list_i = r(r>c);
        list_j = c(r>c);

        particle_i = list_i(~ ismember(list_i,list_j));
        particle_j = list_j(~ ismember(list_i,list_j));

        unique_part = unique(particle_i);

        new_part_centroid = zeros(length(unique_part),2);
        new_part_area = zeros(length(unique_part),1);
        for j0 = 1:length(unique_part)
            part = unique_part(j0);
            associated_part = particle_j(particle_i == part);
            list_particle = cat(1,part,associated_part);

            % computes mass center coordinates 
            new_part_centroid(j0,:) = sum(stats.Area(list_particle).*stats.Centroid(list_particle,:))...
                /sum(stats.Area(list_particle));
            new_part_area(j0) = sum(stats.Area(list_particle));
        end 
        % Get particles that do not need to be merged

        idx_list = (1:size(D,1));
        non_merged_mask = ~ (ismember(idx_list,particle_i) + ismember(idx_list,particle_j));
        non_merged = idx_list(non_merged_mask);

        non_merged_centroid = stats.Centroid(non_merged,:);
        non_merged_area = stats.Area(non_merged);
        new_centroid = cat(1,non_merged_centroid,new_part_centroid);
        new_area = cat(1,non_merged_area,new_part_area);
        
        Centroid = new_centroid;
        Area = new_area;
        new_table = table(Centroid,Area);
        
    else 
        
        new_table = stats;
    end 
    
    new_table = sortrows(new_table,'Centroid','ascend');
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    
    %
%     for i0 = 1:length(list_i)
%         idx_i = list_i(i0);
%         idx_j = list_j(i0);
% 
%         new_table(i0,:) = (stats.Centroid(idx_i,:) + stats.Centroid(idx_j,:))/2;
%     end 
%     
%     mask = stats.Area > min_area;
%     Nb_objects = numel(stats.Area(mask));
%     disp(['Number of detected objects,after erosion :' num2str(Nb_objects)])
%     
    % get position scaled in meter 
%     X_buoys = zeros(Nb_objects,1);
%     Y_buoys = zeros(Nb_objects,1);
    t = (i - 1)/fps; % time 
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
save(file_param,'R0','B0','G0','threshold','SE','min_area','min_distance','-v7.3')

% #############################
%% Extraction for mesange 
% #############################

% Load the whole video 
base = 'E:/Rimouski_2024/Data/2024/0211/Drones/mesange/stereo_001/';
% 'E:/Rimouski_2024/Data/2024/0211/Drones/Fulmar/FULMAR_vertical/';
% addpath(base);

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
f = 2830; % focale in pixel 
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
% create a csv file and write heads 
csv_file = [base 'Buoys_tracking.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,'];
writematrix(heads,csv_file,'WriteMode','append')

% create a csv file for positions in pixels 
csv_file_pix = [base 'Buoys_tracking_pix.csv'];
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

%     img = readFrame(v);

    % Build a distance from colours

    img_double = im2double(img);
    % build intensity matrix 
    I = sqrt((img_double(:,:,1) - R0).^2 + (img_double(:,:,2) - G0).^2 + (img_double(:,:,3) - B0).^2);

    % Mask large values 
    Inew = ones(size(I)) - I;

    Inew(Inew >= threshold) = 1;
    Inew(Inew < threshold) = 0;

    % Detect ROI and extract there position 
    
%     CC = bwconncomp(Inew);
%     disp(['Number of detected objects :' num2str(CC.NumObjects)])
    Inew = imerode(Inew,SE); % erosion 
    CC = bwconncomp(Inew);

    stats = regionprops("table",CC);
% 
%     figure, 
%     imagesc(Inew)
%     hold on 
%     plot(stats.Centroid(:,1),stats.Centroid(:,2),'or')
%     disp(stats)
%  
%     stats.Centroid(:,1) = stats.Centroid(:,1);
%     stats.Centroid(:,2) = stats.Centroid(:,2);
    stats = sortrows(stats,'Centroid','descend'); % sort rows according to x position
    
    
    % Compute distance matrix between objects 
    D = squareform(pdist(stats.Centroid)); 
   
    % Merge particles that are too close to each other 
    mask_dist = (D < min_distance) & (D > 0);
    
    if ismember(1,mask_dist)
        [r,c] = find(mask_dist); % find indices of non-zero elements
        % select only superior part of the matrix 
        list_i = r(r>c);
        list_j = c(r>c);

        particle_i = list_i(~ ismember(list_i,list_j));
        particle_j = list_j(~ ismember(list_i,list_j));

        unique_part = unique(particle_i);

        new_part_centroid = zeros(length(unique_part),2);
        new_part_area = zeros(length(unique_part),1);
        for j0 = 1:length(unique_part)
            part = unique_part(j0);
            associated_part = particle_j(particle_i == part);
            list_particle = cat(1,part,associated_part);

            % computes mass center coordinates 
            new_part_centroid(j0,:) = sum(stats.Area(list_particle).*stats.Centroid(list_particle,:))...
                /sum(stats.Area(list_particle));
            new_part_area(j0) = sum(stats.Area(list_particle));
        end 
        % Get particles that do not need to be merged

        idx_list = (1:size(D,1));
        non_merged_mask = ~ (ismember(idx_list,particle_i) + ismember(idx_list,particle_j));
        non_merged = idx_list(non_merged_mask);

        non_merged_centroid = stats.Centroid(non_merged,:);
        non_merged_area = stats.Area(non_merged);
        new_centroid = cat(1,non_merged_centroid,new_part_centroid);
        new_area = cat(1,non_merged_area,new_part_area);
        
        Centroid = new_centroid;
        Area = new_area;
        new_table = table(Centroid,Area);
        
    else 
        
        new_table = stats;
    end 
    
    new_table = sortrows(new_table,'Centroid','descend');
    Nb_objects = length(new_table.Area(new_table.Area > min_area));
    disp(['Detected objects ' num2str(Nb_objects)])
    
    %
%     for i0 = 1:length(list_i)
%         idx_i = list_i(i0);
%         idx_j = list_j(i0);
% 
%         new_table(i0,:) = (stats.Centroid(idx_i,:) + stats.Centroid(idx_j,:))/2;
%     end 
%     
%     mask = stats.Area > min_area;
%     Nb_objects = numel(stats.Area(mask));
%     disp(['Number of detected objects,after erosion :' num2str(Nb_objects)])
%     
    % get position scaled in meter 
%     X_buoys = zeros(Nb_objects,1);
%     Y_buoys = zeros(Nb_objects,1);
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
file_param = [base 'Parameters_buoy_tracking_mesange_0211'];
save(file_param,'R0','B0','G0','threshold','SE','min_area','min_distance','-v7.3')

% ########################################
%% METHODOLOGY TO EXTRACT COLOUR OF A BUOY 
% ########################################

% fig = figure; 
% imshow(img)
% 
% % coordinates of buoys
% [xi,yi] = getpts(fig);

%%
% % select coordinates of a single buoy
% x0 = round(xi(2));
% y0 = round(yi(2));
% w = 50;
% crop_img = img(y0-w : y0 +w,x0 - w : x0 + w ,:);
% 
% % select pixels of buoy
% % clear xb yb
% % fig_buoy = figure; 
% % imshow(crop_img);
% % [xb,yb] = getpts(fig_buoy);
% 
% P = impixel(im2double(crop_img));
% %%
% % Average of the colours
% idx_x = round(xb);
% idx_y = round(yb);
% 
% c = zeros(3,length(idx_x));
% for i0 = 1:length(idx_x)
%     c(:,i0) = squeeze(crop_img(idx_y(i0),idx_x(i0),:)); % colours RGB
% end 
% 
% mean_RGB = mean(c,2);
% disp(mean_RGB)