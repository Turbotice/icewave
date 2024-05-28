
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
min_area = 5;

% create a csv file and write heads 
csv_file = [base 'Buoys_tracking.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,','X6,','Y6,'];
writematrix(heads,csv_file,'WriteMode','append')

% create a csv file for positions in pixels 
csv_file_pix = [base 'Buoys_tracking_pix.csv'];
heads = ['idx,','t,','x1,','y1,','x2,','y2,','x3,','y3,','x4,','y4,','x5,','y5,','x6,','y6,'];
writematrix(heads,csv_file_pix,'WriteMode','append')
%%

for i = 4500:v.NumFrames

    disp(i)
    % Extract frame i from the video 
    img = read(v,i);
%     img = readFrame(v);

    % Build a distance from colours

    img_double = im2double(img);
    % build intensity matrix 
    I = sqrt((img_double(:,:,1) - R0).^2 + (img_double(:,:,2) - G0).^2 + (img_double(:,:,3) - B0).^2);

    % figure, 
    % imagesc(I)
    % colorbar()
    % figure, 
    % imshow(crop_img),
    % axis image

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
%     figure, 
%     imagesc(Inew)
%     hold on 
%     plot(stats.Centroid(:,1),stats.Centroid(:,2),'or')
%     disp(stats)
 
%     stats.Centroid(:,1) = stats.Centroid(:,1);
%     stats.Centroid(:,2) = stats.Centroid(:,2);
    stats = sortrows(stats,'Centroid','descend'); % sort rows according to x position
    
    
    % Compute distance matrix between objects 
    D = squareform(pdist(stats.Centroid)); 
    min_distance = 20; % minimal distance in pixels to merge objects
    
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
        for i0 = 1:length(unique_part)
            part = unique_part(i0);
            associated_part = particle_j(particle_i == part);
            list_particle = cat(1,part,associated_part);

            % computes mass center coordinates 
            new_part_centroid(i0,:) = sum(stats.Area(list_particle).*stats.Centroid(list_particle,:))...
                /sum(stats.Area(list_particle));
        end 
        % Get particles that do not need to be merged

        idx_list = (1:size(D,1));
        non_merged_mask = ~ (ismember(idx_list,particle_i) + ismember(idx_list,particle_j));
        non_merged = idx_list(non_merged_mask);

        non_merged_centroid = stats.Centroid(non_merged,:);
        new_centroid = cat(1,non_merged_centroid,new_part_centroid);
        new_centroid = sortrows(new_centroid,1,'descend');
        
    else 
        
        new_centroid = stats.Centroid;
    end 
    
    Nb_objects = size(new_centroid,1);
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
    for j = 1:Nb_objects
        xc = new_centroid(j,1);
        yc = new_centroid(j,2);
        % convert coordinates in meters 
        X_buoy = xc/fx_pix;
        Y_buoy = yc/fx_pix; 
        row = cat(2,row,X_buoy,Y_buoy);
        row_pix = cat(2,row_pix,xc,yc);
    end 
    
    row = string(row);
    row_pix = string(row_pix);
    % Write positions of detected ROIs in a csv file 
    writematrix(row,csv_file,'WriteMode','append')
    writematrix(row_pix,csv_file_pix,'WriteMode','append')

end 
disp('Done.')

%% Extraction for bernache 

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
min_area = 5;

% create a csv file and write heads 
csv_file = [base 'Buoys_tracking.csv'];
heads = ['idx,','t,','X1,','Y1,','X2,','Y2,','X3,','Y3,','X4,','Y4,','X5,','Y5,','X6,','Y6,'];
writematrix(heads,csv_file,'WriteMode','append')

% create a csv file for positions in pixels 
csv_file_pix = [base 'Buoys_tracking_pix.csv'];
heads = ['idx,','t,','x1,','y1,','x2,','y2,','x3,','y3,','x4,','y4,','x5,','y5,','x6,','y6,'];
writematrix(heads,csv_file_pix,'WriteMode','append')

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