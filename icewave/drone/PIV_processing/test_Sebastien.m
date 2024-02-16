%%
%##########################################################################
% This code is used to postprocess a data file from a PIV analysis computed
% thanks to PIVlab. It also enables us to generate a structure from the
% postprocessed datas and save it under a .mat file
%##########################################################################

%% Post-process raw datas from PIV analysis and create an associated structure 

% ####### ATTENTION ########
% Be careful to modify values of s and p 
% ##########################

% % Standard PIV Settings
% s = cell(11,2); % To make it more readable, let's create a "settings table"
% %Parameter                          %Setting           %Options
% s{1,1}= 'Int. area 1';              s{1,2}=256;         % window size of first pass
% s{2,1}= 'Step size 1';              s{2,2}=64;         % step of first pass
% s{3,1}= 'Subpix. finder';           s{3,2}=2;          % 1 = 3point Gauss, 2 = 2D Gauss
% s{4,1}= 'Mask';                     s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
% s{5,1}= 'ROI';                      s{5,2}=[3000,1,840,2159];         % Region of interest: [x,y,width,height] in pixels, may be left empty
% %s{5,1}= 'ROI';                      s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
% %s{5,1}= 'ROI';                      s{5,2}=[1,1,895,603];         % Region of interest: [x,y,width,height] in pixels, may be left empty
% s{6,1}= 'Nr. of passes';            s{6,2}=4;          % 1-4 nr. of passes
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
% p{3,1}= 'CLAHE size';            p{3,2}=64;         % CLAHE window size
% p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
% p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
% p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
% p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denoise filter, 0 = disable
% p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size
% p{9,1}= 'Minimum intensity';     p{9,2}=0.0;        % Minimum intensity of input image (0 = no change)
% p{10,1}='Maximum intensity';     p{10,2}=1.0;       % Maximum intensity on input image (1 = no change)


% get the filename 
% date = '20230310';
% base = 'W:/Banquise/Rimouski_2023/Data/drone/';

clear all;
% main directory 
base = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/Traitement_drone_20230310/matData/';
% base = 'E:/PIVlab_drone/matdata/raw_datas/';
folder = [base 'raw_datas_DJI_0308/'];% folder of raw datas
filename = 'PIV_processed_Dt4_b1_W64_full_test_images.mat';
fullname = [folder filename];% filename of raw datas
fx = 1; % spacial sampling frequency (leave it to 1)
ft = 1; % temporal sampling frequency (leave it to 1)
a = 1; % number of boxes to crop on the side
% w = 32; % size of the last window used during PIV process
N = 0;

% directory where we can save the post-processed datas
% save_post_directory = [base 'post_processed/Post_processed_11_02_2023/'];
% if ~exist(save_post_directory)
%     mkdir(save_post_directory)
% end

% load raw_datas
% raw_data_path = [directory filename];
disp(fullname);
disp('Loading Raw Data..');
load(fullname,'u','v','s','p');
disp('Raw Data Loaded');

nb_pass = s{6,2}; % number of pass during the PIV analysis
w = 2^(9-nb_pass); % size of the window during last pass
disp(['Window size = ' num2str(w) ''])
% post-processing
[u_filt,v_filt] = PIV_banquise_postprocessing(u,v,w,N);

%post_pro_fullname = [ base ]

% create a matlab structure from the post-processed data file 
m = genere_structure_banquise(u,v,u_filt,v_filt,fx,ft,a,p,s,w);

% save the structure under a new file name
save_name = filename;
save_name = replace(save_name,'.mat','_total_processed.mat');
savemat_dir = [base 'test/'];
if ~exist(savemat_dir)
    mkdir(savemat_dir)
end
savemat = [savemat_dir save_name];
save(savemat, "m");
disp('DONE.');

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
