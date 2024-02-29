%% Clear 
clear all;
close all;

%% Loading data
date = '20230310';
base = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/Data_ice_Dt4_W64/';
%base = 'W:/Banquise/Rimouski_2023/Data/drone/';
%base = 'E:/Stage MIZ/PIVlab_drone';
%base = '/media/turbots/My Passport//Stage MIZ/PIVlab_drone';

prefixe = 'Stephane_PIV_processed_';
suffixe = 'Dt4_b1_DJI_0402_images_post_processed';
filename = [base prefixe suffixe '.mat'];
%%
disp('Loading Data..');
load(filename);
disp('Data loaded');
%% Creates a folder where to save the generated data
save_folder = base;
if ~exist(save_folder)
    mkdir(save_folder)
end

%%
[ny,nx,nt] = size(m.Vx);

% create a meshgrid
y = (1:1:ny);
x = (nx:-1:1);

[Y,X]=meshgrid(y,x);

% compute the mean component of the velocity field for each frame
Vxmoy = mean(mean(m.Vx,2),1);
Vymoy = mean(mean(m.Vy,2),1);
% reduce the velocity field from its mean value
m.Vx = m.Vx - Vxmoy;
m.Vy = m.Vy - Vymoy;

%% Get rid off quadratic noise (drone movements) on the whole image
[nx,ny,n] = size(m.Vx);
%imax = 105; % maximum index on the x-axis
imax = 105; % maximum idx used for ice 
%imax = nx;
imin = 1; % minimum index on the x-axis
%imin = 25; % minium idx used for waves 
[Y_quad,X_quad] = meshgrid(m.y,m.x);
quadratic_boolean = 1;
% clear Vx_s;

% The quadradtic motion are computed on the whole image !
Vx_s = zeros(nx,ny,n);
for i=1:n
    if quadratic_boolean
        Vx = m.Vx(:,:,i);
        % fit by a quadratic function
        P = polyFit2D(Vx,X_quad,Y_quad,2,2);
        Pth = polyVal2D(P,X_quad,Y_quad,2,2);

        %size(Pth)
        %imagesc(Vx'-Pth')
        Ptot(i,:) = P;
        Vx_s(:,:,i) = Vx-Pth; % get rid off quaratic noise (drone movements)
        %disp('Reduction of quadratic noise');
    else 
        Vx_s(:,:,i) = m.Vx(:,:,i);
    end
end

% We select only the velocity field of the ice 
% V = Vx_s(imin:imax,:,:);

%% Compute and save the smooth temporal FFT 
nb_cuts = 50;
add_pad2 = 2;
disp('Computing FFT_cut')
[TF_cut_t] = compute_time_FFT_cut(Vx_s,nb_cuts,add_pad2);
disp('FFT_cut computed')

% save variables 
matname = [save_folder 'FFT_cut_' num2str(nb_cuts) '_addpad_' num2str(add_pad2) '.mat'];
disp('Saving data..')
save(matname,'TF_cut_t','nb_cuts','add_pad2','-v7.3')
disp('Data saved, work DONE.')