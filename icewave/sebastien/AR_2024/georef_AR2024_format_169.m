%% Script aims at georeferencing images obtained from drone aerial observations

year = '2024';
date = '0921';
drone_ID = 'bernache';
path2data =  ['C:/Users/sebas/Desktop/Amundsen_RA_2024/Data/' year '/' date '/Drones/' drone_ID];

filelist = dir([path2data '/' '*ortho*']);
img_format = '*.JPG';
feet2meter = 1/3.2808; % 1 feet -> x meter

% fields to be converted
keys = {'h','alpha'}; % keys for which unit conversion is needed
initial_units = {'ft','deg'}; % initial units of drone parameters
IS_units = {'m','rad'}; % IS units of drone parameters
coef_conversion = [feet2meter pi/180]; % conversion coefficients from initial unit to IS units


for folder_idx = 1 : 1
    current_folder = filelist(folder_idx).name;
    path_folder = [filelist(folder_idx).folder '\' current_folder];

    % load drone parameters in a structure 
    filelist_drone_param = dir([path_folder '\' '*.txt']);
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
        filelist_image = dir([path_folder '\' img_format]);
        for idx = 1 : length(filelist_image)
            img_name = [filelist_image(idx).folder '\' filelist_image(idx).name];
            disp(img_name)
           
        end 

    end 

end 

%% Try to georeference the image

% define pixel coordinates 
Lx = size(current_img,2);
Ly = size(current_img,1);
% x_pix = (1:1:Lx);
% y_pix = (1:1:Ly);
% [X_pix,Y_pix] = meshgrid(x_pix,y_pix);
% X_pix = permute(X_pix, [2,1]);
% Y_pix = permute(Y_pix, [2,1]);

X_pix=repmat((1:Lx),Ly,1);
Y_pix=repmat((1:Ly)',1,Lx);
% define camera pixel center 
x0 = (Lx + 1)/2;
y0 = (Ly + 1)/2;

% Compute coordinates in real space from pixel coordinates 

%##########################################################################
% I need to compute the focale for format 16:9 -> take a picture of the
% Amundsen and compute the distance -> adjust focale f until it seems ok ? 
% #########################################################################

focale = 1000;
[Xreal,Yreal] = projection_real_space(X_pix,Y_pix,x0,y0,s.h,s.alpha,focale);

figure(1)
imshow(current_img)

img_gray = rgb2gray(current_img);
figure(2),
surf(Xreal,Yreal,img_gray)
shading interp
colormap(gray)
view(2)



%% FUNCTION SECTION 


disp('Function section executed !')