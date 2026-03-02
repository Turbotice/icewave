clear all 
close all

%% Load images 
date = '0226';
drone_ID = 'mesange';
exp_ID = '12-FRAC_001';

base = ['K:/Share_hublot/Data/PIV_images/' date '/Drones/' drone_ID '/' exp_ID '/'];
% base = ['/media/turbots/Elements/Share_hublot/Data/PIV_images/' date '/Drones/' drone_ID '/' exp_ID '/'];
files = dir(fullfile(base,'*.tiff'));

% create cell of full paths
tiffFullPaths = fullfile(base,{files.name});

%% Perform buoys tracking

% define fps
fps = 30;

% define file2save
base_save = ['K:/Share_hublot/Data/' date '/Drones/' drone_ID '/matData/' exp_ID '/'];
file2save = [base_save 'Buoys_tracking_pix_' date '_' drone_ID '_' exp_ID '.mat'];

% set values of RGB channel
R0 = 0.74; % value on canal R
G0 = 0.18; % value on canal G 
B0 = 0.25; % value on canal B 
RGB0 = [R0,G0,B0];

% threshold for binarization 
threshold = 0.78;

% define a structuring element for erosion 
SE = strel('square',2);
min_area = 9; % minimal area of detected objects 
min_distance = 20; % minimal distance in pixels to merge objects

% create an empty structure based on heads
Nb_frames = length(tiffFullPaths);
Nan_array = NaN(Nb_frames,1); % create an array full of Nan values
heads = {'idx','t','x1','y1','x2','y2','x3','y3'};
S_buoys = struct();

% for i = 1:length(heads)
%     S_buoys.(heads{i}) = Nan_array;
% end

% for i0 = 1:Nb_frames
i0 = 7810;
img = imread(tiffFullPaths{i0});
% figure(1);
% imshow(img)

new_table = detect_buoys(img,RGB0,threshold,SE,min_distance);
Nb_objects = length(new_table.Area(new_table.Area > min_area));
disp(['Detected objects ' num2str(Nb_objects)])

t = (i0 - 1)/fps; % time 
S_buoys.idx(i0) = i0 - 1;
S_buoys.t(i0) = t;

% loop over detected objects to fill structure
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
% end

% save(file2save,'S_buoys','-v7.3')
% disp('Done.')