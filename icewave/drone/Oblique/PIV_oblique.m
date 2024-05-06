
%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

base = 'E:/Rimouski_2024/Data/2024/0211/Drones/bernache/matData/stereo_001/';
%base = 'E:/PIVlab_drone/matdata/DJI_0308_Dt4_W64_full/';

filename = 'PIV_processed_N15500_i011500_Dt2_b1_W32_full_post_processed.mat';
matname = [base filename];
%%
disp('Loading Data..');
load(matname);
disp('Data loaded');

%%

Dt = 2; % step between two frames that were compared during the PIV algorithm 
W = 32; 

H = 90;
alpha_0 = 60*pi/180;
fps = 30;
%% Load programm 
base_param = 'C:/Users/sebas/git/icewave/drone/Oblique/';
param_file = [base_param 'Param_Debut_PIV_DJI_20240211153234_0273_D.mat'];

load(param_file)

%%
Vy = m.Vy; 
%%
x_pix=repmat([size(Vy,1):-1:1]*W/2 + 1,size(Vy,2),1)';
y_pix=repmat([size(Vy,2):-1:1]'*W/2 + 1,1,size(Vy,1))';


%%

f=2830;%f=2830; %Focale apparente en pixels

H=param.H; % Hauteur en mètres
alpha_0=param.alpha_0*pi/180; %angle passé en radians
tanalpha_0=tan(alpha_0);
% Les distances sont à prendre par rapport au centre du capteur

x_0=3840/2; % correction à faire, demi-pixel !!
y_0=2160/2;


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = +H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));

tanbeta=(y_pix-y_0)/f;

tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);

X= (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Remise à l'endroit

X=-X;
Y=-Y;

%% 
    figure, 
for i = 1:1: size(Vy,3)
    y_prime = y_pix + squeeze(Vy(:,:,i));
    %y_prime = y_pix;
    Vz = H*(1 - (tanalpha_0*f - (y_prime - y_0)).*(Y./H + 1./tanalpha_0)./(f + tanalpha_0.*(y_prime - y_0)));

    dz(:,:,i) = Vz*fps/Dt;
% 
%     surf(X,Y,dz(:,:,i))
%     shading interp
%     view(2)
%     colorbar()
%     caxis([-1 1])
%     title(num2str(i))
%     pause(0.01)
    
end


%%
figure, 
plot(m.t,squeeze(dz(118,66,:)))

%% 
figure, 
surf(X,Y,dz(:,:,600))
shading interp
view(2)
colorbar()
caxis([-1 1])
title(num2str(i))