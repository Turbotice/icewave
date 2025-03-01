%%
%##########################################################################
% This code is used to postprocess a data file from a PIV analysis computed
% thanks to PIVlab. It also enables us to generate a structure from the
% postprocessed datas and save it under a .mat file
%##########################################################################

%% Post-process raw datas from PIV analysis and create an associated structure 

clear all;
date = '20240223';
base = ['/media/turbots/Hublot24/Share_hublot/Data/' date(5:end) '/Drones/mesange/'];

% base = 'E:/PIVlab_drone/matdata/raw_datas/';
folder = [base 'matData/frac+waves_002/'];% folder of raw datas
filename = 'PIV_processed_i0150_Dt4_b1_W32.mat';
fullname = [folder filename];%[folder filename];% filename of raw datas

a = 1; % number of boxes to crop on the side
w = 32; % size of the last window used during PIV process
Dt = 4; % step between two frames that were compared during the PIV algorithm 
N = 0; % total number of frames processed
i0 = 150; % first index of the frame processed
b = 1; % step between frame A and A' at which velocity is computed

% ##################################################################
%% ################## Create scales ################################
% ##################################################################

facq_t = 29.97; % Frame rate in Hz
ft = 1/facq_t ; % factor scaling for time in sec / frame

% ##########################################
L_x = 3840; % size of the image in pixel, along larger axis (x-axis)
h_drone = 8.53; % height of the drone in meter
theta_x = 34.15; % semi AFOV of the drone, along x-axis, in °
alpha_0 = 90;
facq_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
facq_x = facq_pix*2/w; % scale in box / meter
fx = 1/facq_x; % factor scaling in meter / box
scale_V = (facq_t/Dt) / facq_x; % scale of the velocity in m/s
% ##########################################

% Create t0_UTC
Y = 2024;
M = 02;
D = 11;
H = 15;
MIN = 35;
S = 11;
MS = 690;
TimeZone = 'America/Montreal';
% initial time of recording
t0_UTC = datetime(Y,M,D,H,MIN,S,MS,'TimeZone',TimeZone); 
t0_UTC.TimeZone = 'UTC'; % converts time to UTC time 
t0_UTC.Format = 'yyyy-MMM-dd HH:mm:ss.SSS';

%% load raw_datas

% raw_data_path = [directory filename];
disp(fullname);
disp('Loading Raw Data..');
load(fullname,'u','v','s','p');
disp('Raw Data Loaded');

disp(['Window size = ' num2str(w) ''])
% post-processing
[u_filt,v_filt] = PIV_banquise_postprocessing(u,v,w,N);

%post_pro_fullname = [ base ]

% create a matlab structure from the post-processed data file 
m = genere_structure_banquise(u,v,u_filt,v_filt,fx,ft,scale_V,a,p,s,w,fullname);

% Store parameters used for PIV processing
m.Dt = Dt;
m.i0 = i0;
m.b = b;
m.h_drone = h_drone;
m.theta_x = theta_x;
m.alpha_0 = alpha_0; % pitch angle of the camera in degrees
m.facq_pix = facq_pix;
m.t0_UTC = t0_UTC;

m.units.h_drone = 'm';
m.units.theta_x = 'deg';
m.units.alpha_0 = 'deg';
m.units.facq_pix = 'pix/m';

% save the structure under a new file name
save_name = fullname;
save_name = replace(save_name,'.mat','_total_processed.mat');
savemat_dir = folder;
if ~exist(savemat_dir)
    mkdir(savemat_dir)
end
save([savemat_dir save_name],'m','-v7.3');
disp('Unscaled structure saved');

% scale data 
m = scaling_structure(m); 
disp('Data scaled')

save_name = fullname;
save_name = replace(save_name,'.mat','_scaled.mat');
savemat_dir = folder;
if ~exist(savemat_dir)
    mkdir(savemat_dir)
end
save([savemat_dir save_name],'m','-v7.3');
disp('Scaled structure saved');

base = 'W:/SagWin2024/Data/0223/Drones/bernache/matData/12-waves_010/';
filename = [base 'physical_parameters.txt'];
parameters_list = {'h_drone','theta_x','alpha_0','facq_pix','fx','ft','scale_V'};

saving_parameters(m,filename,parameters_list)
disp('Parameters saved in a txt file')

disp('DONE.')
%%
% ###############################################
% %% Try parameters of the PIV between two pictures 
% % ###############################################
% 
% % Standard PIV Settings
% s = cell(11,2); % To make it more readable, let's create a "settings table"
% %Parameter                          %Setting           %Options
% s{1,1}= 'Int. area 1';              s{1,2}=128;         % window size of first pass
% s{2,1}= 'Step size 1';              s{2,2}=32;         % step of first pass
% s{3,1}= 'Subpix. finder';           s{3,2}=2;          % 1 = 3point Gauss, 2 = 2D Gauss
% s{4,1}= 'Mask';                     s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
% %s{5,1}= 'ROI';                      s{5,2}=[1,1080,1023,500];         % Region of interest: [x,y,width,height] in pixels, may be left empty
% s{5,1}= 'ROI';                      s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
% %s{5,1}= 'ROI';                      s{5,2}=[1,1,895,603];         % Region of interest: [x,y,width,height] in pixels, may be left empty
% s{6,1}= 'Nr. of passes';            s{6,2}=3;          % 1-4 nr. of passes
% s{7,1}= 'Int. area 2';              s{7,2}=128;         % second pass window size
% s{8,1}= 'Int. area 3';              s{8,2}=64;         % third pass window size
% s{9,1}= 'Int. area 4';              s{9,2}=32;         % fourth pass window size
% s{10,1}='Window deformation';       s{10,2}='*spline'; % '*spline' is more accurate, but slower
% s{11,1}='Repeated Correlation';     s{11,2}=0;         % 0 or 1 : Repeat the correlation four times and multiply the correlation matrices.
% s{12,1}='Disable Autocorrelation';  s{12,2}=0;         % 0 or 1 : Disable Autocorrelation in the first pass.
% s{13,1}='Correlation style';        s{13,2}=0;         % 0 or 1 : Use circular correlation (0) or linear correlation (1).
% s{14,1}='Repeat last pass';   s{14,2}=0; % 0 or 1 : Repeat the last pass of a multipass analyis
% s{15,1}='Last pass quality slope';   s{15,2}=0.025; % Repetitions of last pass will stop when the average difference to the previous pass is less than this number.
% 
% % Standard image preprocessing settings
% p = cell(10,1);
% %Parameter                       %Setting           %Options
% p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
% p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
% p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
% p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
% p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
% p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
% p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denoise filter, 0 = disable
% p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size
% p{9,1}= 'Minimum intensity';     p{9,2}=0.0;        % Minimum intensity of input image (0 = no change)
% p{10,1}='Maximum intensity';     p{10,2}=1.0;       % Maximum intensity on input image (1 = no change)
% 
% 
% %% Check if the Dt is coherent for PIV 
% 
% date = '20230310';
% base = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Data/drone/';
% 
% % arguments of piv_analysis
% directory = [base date '/contexte/video/DJI_0402_images/'];
% filename1 = 'im_0000.tiff';
% filename2 = 'im_0003.tiff';
% nr_of_cores = 20;
% graph = 0;
% 
% [x, y, u, v, typevec,corr_map] = piv_analysis(directory, filename1, filename2,p,s,nr_of_cores,graph);
% 
% %% Plot velocity field 
% 
% [nx,ny] = size(u);
% 
% % create a meshgrid
% ymesh = (1:1:ny);
% %xmesh = (nx:-1:1);
% xmesh = (1:1:nx);
% 
% [Y,X]=meshgrid(ymesh,xmesh);
% % compute the mean component of the velocity field
% %{Vxmoy = mean(mean(m.Vx,2),1);
% Vymoy = mean(mean(m.Vy,2),1);
% % reduce the velocity field from its mean value
% m.Vx = m.Vx - Vxmoy;
% m.Vy = m.Vy - Vymoy;
% %}
% 
% figure(1)
% 
% subplot(2,1,1)
% hold off
% surf(Y,X,u)
% shading interp
% view(2)
% caxis([-10 10])
% colorbar()
% title('Vx')
% xlabel('y'), ylabel('x')
% 
% subplot(2,1,2)
% hold off
% surf(Y,X,v)
% 
% shading interp
% view(2)
% caxis([-10 10])
% colorbar()
% title('Vy')
% xlabel('y'), ylabel('x')
% 
% 
% % try a quiver plot
% figure(2)
% 
% quiver(Y,X,u,v)
% 
% %% Show an image and try to select the correct zone 
% 
% date = '20230310';
% base = 'W:/Banquise/Rimouski_2023/Data/drone/';
% 
% % arguments of piv_analysis
% directory = [base date '/contexte/video/DJI_0402_images/'];
% filename = 'im_0000.tiff';
% 
% file = [directory filename];
% image = imread(file);
% [Ny ,Nx , Nc] = size(image);
% yroi = 1;
% heightroi = Ny-1;
% xroi = 3000;
% widthroi = Nx-1 - xroi;
% new_image = image(yroi:yroi+heightroi,xroi:xroi+widthroi);
% 
% figure(3)
% imshow(image);
% 
% figure(4)
% 
% imshow(new_image)
% 
