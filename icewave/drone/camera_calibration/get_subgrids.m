% This script enables detection of subgrids within a picture of tiling, 
% checkerboard or other regular square pattern
% 
% Author : Sebastien Kuchly, PMMH, ESPCI Paris

% 1) we load binarized images, obtained from im_2_binary.m
% 2) squares are detected on the binarized image, function
% calib_drone_getsquares
% 3) manual detection of subgrids, a file containing all subgrids from a
% single image is saved, function calib_drone_getsubgrid
% 4) all subgrids are saved in a single matrix
%

clear all
close all

%% Load binarised images and build main structure

M={};

base = 'W:';
main_path = [ base '/Banquise/Calibration_PMMH_2026/20251216/calib_video_4K/frame_selection/'];
folder = [main_path 'Binarisation/'];

% collect list of images and binarized images 
bwlist = dir([folder '*.mat']);
imlist = dir([folder '*.tiff']);

figure(1)
for i=1:length(bwlist)
    name = bwlist(i).name;
    imname = imlist(i).name;

    disp(name)
    disp(imname)
    load(fullfile(bwlist(i).folder,name)) % load binarized image
    im = imread(fullfile(imlist(i).folder,imname)); % load image 

    m.name = name;
    m.imname = imlist(i).name;
    M(i).m=m;

end

%% Get tiles of floor tiling

for i=1:length(M)
    name = M(i).m.name;
    imname = M(i).m.imname;
    load(fullfile(folder,name))
    im = imread(fullfile(imlist(i).folder,imname));

    figure(5)
    imshow(BW)
    hold on

    [X,Y] = calib_drone_getsquares(im,BW);
    %    [X,Y]=calib_drone_getgrid(im,BW);

    M(i).m.X=X;
    M(i).m.Y=Y;
    title(name)
    pause(2)
end

%% Collect subgrid of tiles 

save_folder = [main_path 'Transient/'];
if ~isfolder(save_folder)
    disp([save_folder 'has been created'])
    mkdir(save_folder)
end

for i=8:length(M)
    name = M(i).m.name;
    imname = M(i).m.imname;
    load(fullfile(folder,name))
    im = imread(fullfile(imlist(i).folder,imname));

    figure(5)
    imshow(BW)
    hold on

    X = M(i).m.X;
    Y = M(i).m.Y;

    Msub = calib_drone_getsubgrid(im,X,Y,replace(imname,'.','_'),save_folder);
    %    [X,Y]=calib_drone_getgrid(im,BW);

    title(name)
    pause(2)
end

%% Load transient file : contains several subgrids for each image 

Mp = {};

matlist = dir([save_folder 'transient_im_*.mat']);
nx = 9;
ny = 7;
c=0; % counter
for i=1:length(matlist)
    Msub = load([save_folder matlist(i).name]);

    for k=1:length(Msub.M)
        n = length(Msub.M(k).m.Xsub);
        if n == nx*ny
            c=c+1;
            Mp(c).m=Mref(i).m;
            Mp(c).m.X = Msub.M(k).m.Xsub;
            Mp(c).m.Y = Msub.M(k).m.Ysub;
            % Mp(c).m.corners = Msub.M(k).m.corners;
        end
    end
    disp(c)
end

%% Save structure Mp

file2save = [main_path 'matrix_points.mat'];
save(file2save,'Mp','-v7.3')