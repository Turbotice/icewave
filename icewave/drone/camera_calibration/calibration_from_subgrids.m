%% This scrip enables to correctly oriente subgrids in order to use them for
% camera calibration 

% Author : Sebastien Kuchly, PMMH, ESPCI Paris

% 1) we load subgrid points obtained from get_subgrids.m
% 2) we look for each four corners of each subgrid, function find_corners
% 3) points are sorted within each subgrid according to generatecheckerboardpoints ordering
% 4) camera parameters are then estimated
%

%% Load data 

% path2data = ['Z:/skuchly/Codes/' 'matrix_points_bernache.mat'];
main_path = 'W:/Banquise/Calibration_PMMH_2026/20251216/calib_video_4K/frame_selection/';
path2data = [main_path 'matrix_points.mat'];
load(path2data)

%% Show a given subgrid
nx = 9; 
ny = 7;

Mp = find_corners(Mp,nx,ny);

%% Sort all points of a subgrid according to generatecheckerboardpoints ordering

Mtab = sort_points(Mp,nx,ny);

%% Save Mtab

filename = [main_path 'matrix_points_sorted.mat'];

save(filename,'Mtab','-v7.3')

%% Estimate camera calibration points 

% collect imagepoints
for i = 1:length(Mtab)
    imagePoints(:,:,i) = Mtab(i).points;
end

% generate checkerboard pattern 
boardSize = [8,10]; % adapt to size of the detected subgrid
squareSizeInMM = 3955 / 40; % adapt scale
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
imageSize = [3840,2160]; % adapt image size 

% calibration
stereoParams = estimateCameraParameters(imagePoints,worldPoints,ImageSize=imageSize);


%% Save stereoParams and show extrinsics parameters

figure(5)
showExtrinsics(stereoParams)

file_stereo = [main_path 'stereoParams_Tourterelle_4K.mat'];
save(file_stereo,'stereoParams','-v7.3')












%% 

i = 31;
corners = Mp(i).m.corners;
X = Mp(i).m.X;
Y = Mp(i).m.Y;
points = [X(:),Y(:)];
ordered_pts = points;
T = 10;

% find line between corner 1 and 2
[xtop,ytop] = find_line(X,Y,corners(1,1),corners(1,2),...
    corners(2,1),corners(2,2),T);

% find line between corner 4 and 3
[xbottom,ybottom] = find_line(X,Y,corners(4,1),corners(4,2),...
    corners(3,1),corners(3,2),T);


figure(4),
plot(points(:,1),points(:,2),'+')
hold on 
plot(corners(:,1),corners(:,2),'ks')
hold on 
plot(xtop,ytop,'ro')
hold on 
plot(xbottom,ybottom,'dg')




