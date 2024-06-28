
clear all
close all

%% This script details method to georeference an image from drone view 

base = 'W:/SagWin2024/Data/0211/Drones/Initial_images/';
filename = [base 'bernache_im_11500.tiff'];

image_bernache = imread(filename); 

% load parameters 
base_param = 'W:/SagWin2024/Data/0211/Drones/bernache/matData/18-stereo_001/';
param_file = [base_param 'Param_Debut_bernache_PIV_DJI_20240211153234_0273_D.mat'];

load(param_file)

% load objects position 
obj_file = [base_param 'POS_DJI_20240211153234_0273_D.mat'];
load(obj_file)
%% 
H = param.H;
alpha_0 = param.alpha_0*pi/180;
f = 2830; % focale apparente en pixel 

%%
% create array of pixels
x_pix=repmat([1:size(image_bernache,2)],size(image_bernache,1),1);
y_pix=repmat([1:size(image_bernache,1)]',1,size(image_bernache,2));

% Les distances sont à prendre par rapport au centre du capteur

x_0=3840/2;
y_0=2160/2;


L_0 = H/sin(alpha_0);


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = (y_pix - y_0)*H/sin(alpha_0)./(f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0));
X = (x_pix - x_0)*H./(f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0));

% Formules calculées à Rimouski... Fausses ?
% Y = H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));
% 
% tanbeta=(y_pix-y_0)/f;
% tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);
% 
% X= (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Remise à l'endroit

X=X;
Y=-Y;

%%
figure,
surf(X,Y,image_bernache(:,:,1))
Z1=max(image_bernache(:,:,1),[],'all');

colormap(gray)
shading interp
axis equal
view(2)
title('bernache')
% hold on

%%

% Puis les points cliqués

x_points_pix = POS(:,2); %direction horizontale du capteur
y_points_pix = POS(:,3); % direction verticale du capteur


x_0=3840/2*ones(size(x_points_pix));
y_0=2160/2*ones(size(x_points_pix));

 
 % On calcule les distances dans le plan à partir des coordonnées capteurs.
% 

Y_points= L_0 *(y_points_pix-y_0)./(f*sin(alpha_0) + (y_points_pix-y_0) *cos(alpha_0));
X_points = (x_points_pix-x_0)./f.*(L_0 - (cos(alpha_0) *L_0 *(y_points_pix-y_0))./(f*sin(alpha_0) + (y_points_pix-y_0)*cos(alpha_0)));


% % Remise à l'endroit
% 
 X_points = X_points;
 Y_points = -Y_points;
% 
% 
% % % Remise en 0,0 de la caisse noire. Sur l'image, sa position est en X_points(1),Y_points(1) 
% % % On fait
%  X=X-X_points(1);
%  Y=Y-Y_points(1);
% %   
% % % Recalage des points à l'origine de la caisse noire
%  X_points=X_points-X_points(1);
%  Y_points=Y_points-Y_points(1);

% Figure de l'image redressée


figure,
surf(X,Y,image_bernache(:,:,1))
Z1=max(image_bernache(:,:,1),[],'all');

colormap(gray)
shading interp
axis equal
view(2)
title('bernache')
hold on

% 
% % On les ajoute sur le plot
% 
 plot3(X_points,Y_points,1.1*double(Z1)*ones(size(X_points)),'rd','MarkerSize',4,'MarkerFaceColor','r');

 % #################################
 %% ##### Rectification for mesange ######
 % #################################
 
 base = 'W:/SagWin2024/Data/0211/Drones/Initial_images/';
filename = [base 'mesange_im_11496.tiff'];

image_bernache = imread(filename); 

% load parameters 
base_param = 'W:/SagWin2024/Data/0211/Drones/mesange/matData/2-stereo_001/';
param_file = [base_param 'Param_Debut_mesange_PIV_DJI_20240211213233_0222_D.mat'];

load(param_file)

% load objects position 
obj_file = [base_param 'POS_Debut_PIV_DJI_20240211213233_0222_D.mat'];
load(obj_file)
%% 
H = param.H;
alpha_0 = param.alpha_0*pi/180;
f = 2830; % focale apparente en pixel 

%%
% create array of pixels
x_pix=repmat([1:size(image_bernache,2)],size(image_bernache,1),1);
y_pix=repmat([1:size(image_bernache,1)]',1,size(image_bernache,2));

% Les distances sont à prendre par rapport au centre du capteur

x_0=3840/2;
y_0=2160/2;


L_0 = H/sin(alpha_0);


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = (y_pix - y_0)*H/sin(alpha_0)./(f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0));
X = (x_pix - x_0)*H./(f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0));

% Formules calculées à Rimouski... Fausses ?
% Y = H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));
% 
% tanbeta=(y_pix-y_0)/f;
% tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);
% 
% X= (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Remise à l'endroit

X=X;
Y=-Y;

%%
figure,
surf(X,Y,image_bernache(:,:,1))
Z1=max(image_bernache(:,:,1),[],'all');

colormap(gray)
shading interp
axis equal
view(2)
title('bernache')
% hold on

%%

% Puis les points cliqués

x_points_pix = POS(:,2); %direction horizontale du capteur
y_points_pix = POS(:,3); % direction verticale du capteur


x_0=3840/2*ones(size(x_points_pix));
y_0=2160/2*ones(size(x_points_pix));

 
 % On calcule les distances dans le plan à partir des coordonnées capteurs.
% 

Y_points= L_0 *(y_points_pix-y_0)./(f*sin(alpha_0) + (y_points_pix-y_0) *cos(alpha_0));
X_points = (x_points_pix-x_0)./f.*(L_0 - (cos(alpha_0) *L_0 *(y_points_pix-y_0))./(f*sin(alpha_0) + (y_points_pix-y_0)*cos(alpha_0)));


% % Remise à l'endroit
% 
 X_points = X_points;
 Y_points = -Y_points;
% 
% 
% % % Remise en 0,0 de la caisse noire. Sur l'image, sa position est en X_points(1),Y_points(1) 
% % % On fait
%  X=X-X_points(1);
%  Y=Y-Y_points(1);
% %   
% % % Recalage des points à l'origine de la caisse noire
%  X_points=X_points-X_points(1);
%  Y_points=Y_points-Y_points(1);

% Figure de l'image redressée


figure,
surf(X,Y,image_bernache(:,:,1))
Z1=max(image_bernache(:,:,1),[],'all');

colormap(gray)
shading interp
axis equal
view(2)
title('mesange')
hold on

% 
% % On les ajoute sur le plot
% 
 plot3(X_points,Y_points,1.1*double(Z1)*ones(size(X_points)),'rd','MarkerSize',4,'MarkerFaceColor','r');

