%% Script aims at georeferencing images extracted from a drone video of format 4K
package_path = 'C:/Users/sebas/git/icewave/drone/Function_package/';
addpath(package_path, '-end');

year = '2024';
date = '0921';
drone_ID = 'bernache';
exp_ID = '06-waves_001';
path2data =  ['F:/Amundsen_RA_2024/Data/' date '/Drones/' drone_ID '/' exp_ID '/'];

filelist = dir(path2data);
img_format = '*.tiff';
feet2meter = 1/3.2808; % 1 feet -> x meter

% fields to be converted
keys = {'h','alpha'}; % keys for which unit conversion is needed
initial_units = {'ft','deg'}; % initial units of drone parameters
IS_units = {'m','rad'}; % IS units of drone parameters
coef_conversion = [feet2meter pi/180]; % conversion coefficients from initial unit to IS units

% Create a folder where results can be saved
fig_folder = [filelist(1).folder '\' filelist(1).name '\Figures\'];
if ~isfolder(fig_folder)
    mkdir(fig_folder)
end
fig_resolution = 300; % resolution of saved fig

% Parameters for image binarization
default_sensitivity = 0.58; % sensitivity used by adaptthresh 
radius_open = 3; % radius of structuring element to perform erosion
minArea = 100; % minimal area (in pixels) of detected objects
conn = 4; % connection used for bwboundaries

%% Perform ice floes detection
path_folder = path2data;

% load drone parameters in a structure 
filelist_drone_param = dir([path_folder '\' '*parameters.txt']);
for file_idx = 1 : length(filelist_drone_param)
    current_filename = filelist_drone_param(file_idx).name;
    lines = readcell([filelist_drone_param(file_idx).folder '\' current_filename]);
    s = struct();
    s.units = struct();
    for variable_idx = 1: size(lines,1)
        header = lines{variable_idx,1};
        s.(header) = lines{variable_idx,3};
        s.units.(header) = lines{variable_idx,4};
    end 

    % Change drone parameters coefficient 
    s = conversion_structure_units(s,keys,initial_units,IS_units,coef_conversion);

    % load JPG images one by one 
    filelist_image = dir([path_folder img_format]);
    for idx = 1 : length(filelist_image)
        img_name = [filelist_image(idx).folder '\' filelist_image(idx).name];
        disp(img_name)

        current_img = imread(img_name);
        
        txt_cell = split(filelist_image(idx).name,'.');
        prefix_current_fig = txt_cell{1};
       
        BW_param = struct('sensitivity',default_sensitivity,'radius_open',radius_open,...
            'minArea',minArea,'conn',conn);

        results = ice_floes_detection(current_img,s,BW_param,fig_folder,...
            prefix_current_fig,fig_resolution);

        results.figfilename = img_name;

        file2save = [fig_folder prefix_current_fig '_Data_sensitivity' num2str(results.BW_param.sensitivity) ...
            '_imopen' num2str(results.BW_param.radius_open) '_conn' num2str(results.BW_param.conn) '_minarea' ...
            num2str(results.BW_param.minArea) '.mat'];

        disp('Saving Data..')
        save(file2save,'results','-v7.3')
        disp('Data saved !')
       
    end 

end 

