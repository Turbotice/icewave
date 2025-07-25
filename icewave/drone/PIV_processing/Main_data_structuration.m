%%
%##########################################################################
% This code is used to postprocess a data file from a PIV analysis computed
% thanks to PIVlab. It also enables us to generate several files : 
% - a structure from the postprocessed datas and save it under a .mat file
% - a similar structure which contains scaled data
% - a file .txt which contains all relevant parameters 
%##########################################################################

%% Post-process raw datas from PIV analysis and create an associated structure 

clear all;
date = '0226';
drone_name = 'mesange';
exp_ID = '12-FRAC_001';
ID = [date '_' drone_name '_' exp_ID];

base = ['/media/turbots/Elements/Share_hublot/Data/' date '/Drones/' drone_name '/'];
folder = [base 'matData/' exp_ID '/'];% folder of raw datas
filename = 'PIV_processed_i00_N0_Dt7_b1_W32_xROI1_width3839_yROI1_height2159.mat'; % file to load
fullname = [folder filename];

% ##################################################################
%% ################## Input scales and parameters################################
% ##################################################################

% ##########################################
Lx = 3840; % size of the image in pixel, along larger axis (x-axis)
Ly = 2160; % size of the image in pixel, along minor axis (y-axis)
h_drone = 64.6; % height of the drone in meter
alpha_0 = 90*pi/180; % camera pitch angle to the horizontal

x_0 = (Lx + 1)/2; % camera sensor center
y_0 = (Ly + 1)/2; % camera sensor center

focale = 2700; %in pixels 
theta_x = atan(Lx/focale/2); % semi AFOV of the drone, along x-axis, in °
 
facq_pix = Lx/(2*h_drone*tan(theta_x)); % scale in pixels / meter
facq_t = 30; % Frame rate in Hz
ft = 1/facq_t ; % factor scaling for time in sec / frame
% ##########################################

% Longitude and latitude during flight 
latitude = 48.34836;
longitude = -68.81427;

% Create t0 local (beginning of flight)
Y = 2024;
M = 02;
D = 26;
H = 20; % local time of drone
MIN = 22;
S = 00;
MS = 509;
 
if strcmp(drone_name,'mesange') & Y == 2024
    TimeZone = 'Europe/Paris'; % mesange 
else 
    TimeZone = 'America/Montreal'; % bernache or Fulmar
end 

% initial time of recording
t0_UTC = datetime(Y,M,D,H,MIN,S,MS,'TimeZone',TimeZone); 
t0_UTC.TimeZone = 'UTC'; % converts time to UTC time 
t0_UTC.Format = 'yyyy-MM-dd HH:mm:ss.SSS';

% #############################################################
%% ############# Load raw_datas and create structures #########
% #############################################################

% raw_data_path = [directory filename];
disp(fullname);
disp('Loading Raw Data..');
load(fullname,'u','v','x_pix','y_pix','PIV_param');
disp('Raw Data Loaded');

% PIV parameters 
a = 0; % number of boxes to crop on the side
w = PIV_param.w; % size of the last window used during PIV process
Dt = PIV_param.Dt; % step between two frames that were compared during the PIV algorithm 
N = PIV_param.N; % total number of frames processed
i0 = PIV_param.i0; % first index of the frame processed
b = PIV_param.b; % step between frame A and A' at which velocity is computed

disp(['Interrogation window size = ' num2str(w) ''])

% define space scaling for PIV boxes
facq_x = facq_pix*2/w; % scale in box / meter
fx = 1/facq_x; % factor scaling in meter / box
scale_V = (facq_t/Dt) / facq_pix; % scale of the velocity in m.frame/pixel/s

% post-processing
[u_filt,v_filt] = PIV_banquise_postprocessing(u,v,w,N);

% create a matlab structure from the post-processed data file 
m = genere_structure_banquise(u,v,u_filt,v_filt,x_pix,y_pix,a,fullname);

% Create pixel coordinates of image center 
m.PIXEL.x0 = (Lx + 1)/2;
m.PIXEL.y0 = (Ly + 1)/2;

% Store parameters used for PIV processing
m.PIV_param.Dt = Dt;
m.PIV_param.w = w;
m.PIV_param.i0 = i0;
m.PIV_param.N = N;
m.PIV_param.a = a; 
m.PIV_param.b = b;
m.PIV_param.ROI = PIV_param.ROI;
m.PIV_param.p_param = PIV_param.p;
m.PIV_param.s_param = PIV_param.s;

% Store scales
if alpha_0 == pi/2
    m.SCALE.scale_V = scale_V;
    m.SCALE.fx = fx;
    m.SCALE.facq_pix = facq_pix;
    m.SCALE.ft = ft;
    m.SCALE.facq_t = facq_t;

    m.SCALE.units.scale_V = 'm.frame/pixel/s';
    m.SCALE.units.fx = 'm/box';
    m.SCALE.units.facq_pix = 'pix/m';
    m.SCALE.units.ft = 's/frame';
    m.SCALE.units.facq_t = 'frame/s';
else
    m.SCALE.scale_V = [];
    m.SCALE.fx = [];
    m.SCALE.facq_pix = [];
    m.SCALE.ft = ft;
    m.SCALE.facq_t = facq_t;

    m.SCALE.units.scale_V = 'oblique view';
    m.SCALE.units.fx = 'oblique view';
    m.SCALE.units.facq_pix = 'oblique view';
    m.SCALE.units.ft = 's/frame';
    m.SCALE.units.facq_t = 'frame/s';
    
end

% Store Physical parameters 
m.DRONE.h_drone = h_drone;
m.DRONE.theta_x = theta_x;
m.DRONE.alpha_0 = alpha_0; % pitch angle of the camera in degrees
m.DRONE.focale = focale;
m.DRONE.Lx = Lx;
m.DRONE.Ly = Ly;
m.DRONE.name = drone_name;

m.DRONE.units.h_drone = 'm';
m.DRONE.units.theta_x = 'rad';
m.DRONE.units.alpha_0 = 'rad';
m.DRONE.units.focale = 'pix';
m.DRONE.units.Lx = 'pix';
m.DRONE.units.Ly = 'pix';

% Store initial GPS position
m.GPS.latitude = latitude ;
m.GPS.longitude = longitude ;
m.GPS.units.latitude = 'deg';
m.GPS.units.longitude = 'deg';

% store inital UTC_time
m.t0_UTC = t0_UTC;

% store ID experiment
m.ID = ID;

% save the structure under a new file name
save_name = fullname;
save_name = replace(save_name,'.mat','_total_processed.mat');
save(save_name,'m','-v7.3');
disp('Unscaled structure saved');

% scale data 
if alpha_0 == 90*pi/180
    m = scaling_structure(m,m.SCALE.scale_V,m.SCALE.fx,m.SCALE.ft); 
else
    m = scaling_structure_oblique(m,double(m.PIXEL.x_pix),double(m.PIXEL.y_pix),m.t,m.PIXEL.x0,m.PIXEL.y0,...
        m.DRONE.h_drone,m.DRONE.alpha_0,m.DRONE.focale,m.SCALE.facq_t,m.PIV_param.Dt);
end 
% convert UTC time in character
m.t0_UTC = char(t0_UTC);
disp('Data scaled')

save_name = fullname;
save_name = replace(save_name,'.mat','_scaled.mat');
save(save_name,'m','-v7.3');
disp('Scaled structure saved');

filename = [folder 'physical_parameters.txt'];
parameters_list = {'ID'};
saving_parameters(m,filename,parameters_list)
parameters_list = {'name','h_drone','theta_x','alpha_0','focale'}; 
saving_parameters(m.DRONE,filename,parameters_list)
parameters_list = {'scale_V','fx','facq_pix','ft'};
saving_parameters(m.SCALE,filename,parameters_list)
disp('Parameters saved in a txt file')

disp('DONE.')