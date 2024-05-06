clear all; 
close all; 

%% Load PIV 

base = '//192.168.1.70/Share/Data/0211/Drones/bernache/matData/stereo_001/';
%base = 'E:/PIVlab_drone/matdata/DJI_0308_Dt4_W64_full/';

filename = 'PIV_processed_N15500_i011500_Dt2_b1_W32_full_post_processed.mat';
matname = [base filename];

%% Load param inclinaison
base_param = 'C:/Users/sebas/git/icewave/drone/Oblique/';
param_name = 'Param_Debut_PIV_DJI_20240211153234_0273_D.mat';

param_file = [base_param param_name];
load(param_file)


%%
disp('Loading Data..');
load(matname);
disp('Data loaded');

%%
L_x = 3840; % size of the image in pixel, along x-axis
h_drone = 294.6*0.3048; % height of the drone in meter
theta_x = 32.75; % semi AFOV of the drone, along x-axis, in °
alpha_0 = 60.1; 

fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
W = 32; 

%%
Vy = m.Vy;

%% 


x_pix=repmat([size(Vy,1):-1:1]*W/2 + 1,size(Vy,2),1)';
y_pix=repmat([size(Vy,2):-1:1]'*W/2 + 1,1,size(Vy,1))';

%%
figure, 
imagesc(x_pix)
colorbar()


%%


f=2830;%f=2830; %Focale apparente en pixels

H=param.H; % Hauteur en mètres
alpha_0=param.alpha_0*pi/180; %angle passé en radians
tanalpha_0=tan(alpha_0);
% Les distances sont à prendre par rapport au centre du capteur

x_0 = 3840/2; % correction à faire, demi-pixel !!
y_0 = 2160/2;


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));

tanbeta=(y_pix-y_0)/f;

tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);

X = (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Remise à l'endroit

X=-X;
Y=-Y;

%%
figure, 
for i = 1:size(Vy,3)
    
    surf(X,Y,Vy(:,:,i))
    shading interp
    view(2)
    title(num2str(i))
    pause(0.01)

end

%% 

for i0 = 1:5:size(Vy,3)

    y_prime = y_pix + Vy(:,:,i0);

    Vz = H*( 1 - (tanalpha_0*f - (y_prime - y_0))./(f + tanalpha_0*(y_prime - y_0)) .* (Y./H + 1/tanalpha_0));
    dz(:,:,i0) = Vz.*fps/Dt;
    
    disp(i0)
    
    surf(X,Y,squeeze(dz(:,:,i0)))
    shading interp
    colorbar()
    view(2)
    title(num2str(i0))
    pause(0.01)

    
end

%%
x_c = floor(size(dz,2)/2);
y_c = floor(size(dz,2)/2 + 1);

figure, 
plot(m.t,squeeze(dz(x_c,y_c,:)))
