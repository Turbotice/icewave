close all
clear all

%% Chargement des datas

% Path depuis le dossier parent

filepathbernache='bernache/18_stereo_01/Position_ref_ligne';
filepathmesange='mesange/2-stereo_01';
local_path='/Users/eddiaca/Desktop/PMMH/Banquise/Rimouski_2024/Redressement_Saguenay/Test_optimisation_focale/';
% chargement 

cd ..

cd(filepathbernache)
k=rdir('Param*');
g=rdir('POS*.mat');
l=rdir('Debut*.tif');

param_bernache=load(k{1});
POS_objets_bernache=load(g{1});
image_bernache=imread(l{1});

cd(local_path)



cd ..

cd(filepathmesange)

k=rdir('Param*.mat');
g=rdir('POS*.mat');
l=rdir('Debut*.tif');


param_mesange=load(k{1});
POS_objets_mesange=load(g{1});
image_mesange=imread(l{1});

clear g k l

cd(local_path)

%% 

for f=2830%[2830-200:25:2830+200]

%f=2830; %Focale apparente en pixels
close all

%% Redressement bernache

% On commence par l'image

H=param_bernache.param.H; % Hauteur en mètres
alpha_0=param_bernache.param.alpha_0*pi/180; %angle passé en radians
tanalpha_0=tan(alpha_0);


x_pix=repmat([size(image_bernache,2):-1:1],size(image_bernache,1),1);
y_pix=repmat([size(image_bernache,1):-1:1]',1,size(image_bernache,2));

% Les distances sont à prendre par rapport au centre du capteur

x_0=3840/2;
y_0=2160/2;




% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));

tanbeta=(y_pix-y_0)/f;
tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);

X= (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Remise à l'endroit

X=-X;
Y=-Y;



% Puis les points cliqués

x_points_pix=-POS_objets_bernache.POS(:,2); %direction horizontale du capteur
y__points_pix=-POS_objets_bernache.POS(:,3); % direction verticale du capteur




x_0=-3840/2*ones(size(x_points_pix));
y_0=-2160/2*ones(size(x_points_pix));


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y_points = H/tanalpha_0 - H *(tanalpha_0 *(y__points_pix-y_0)+f)./(f*tanalpha_0 - (y__points_pix-y_0));

tanbeta_points=(y__points_pix-y_0)/f;
tanalpha_y_points=(tanalpha_0-tanbeta_points)./(1+ tanalpha_0*tanbeta_points);

X_points= (x_points_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y_points.^2);

% Remise à l'endroit

X_points=-X_points;
Y_points=-Y_points;


% Remise en 0,0 de la caisse noire. Sur l'image, sa position est en X_points(1),Y_points(1) 
% On fait
X=X-X_points(1);
Y=Y-Y_points(1);
  
% Recalage des points à l'origine de la caisse noire
X_points=X_points-X_points(1);
Y_points=Y_points-Y_points(1);

% Figure de l'image redressée


figure(2)
surf(X,Y,image_bernache(:,:,1))
Z1=max(image_bernache(:,:,1),[],'all');

colormap(gray)
shading interp
axis equal
view(2)
title('bernache')
hold on


% On les ajoute sur le plot

plot3(X_points,Y_points,1.1*double(Z1)*ones(size(X_points)),'rd','MarkerSize',18,'MarkerFaceColor','r');


%% Redressement mesange



% On commence par l'image

H=param_mesange.param.H; % Hauteur en mètres
alpha_0=param_mesange.param.alpha_0*pi/180; %angle passé en radians
tanalpha_0=tan(alpha_0);


x_pix=repmat([size(image_mesange,2):-1:1],size(image_mesange,1),1);
y_pix=repmat([size(image_mesange,1):-1:1]',1,size(image_mesange,2));

% Les distances sont à prendre par rapport au centre du capteur

x_0=3840/2;
y_0=2160/2;




% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));

tanbeta=(y_pix-y_0)/f;
tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);

X= (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Pas de remise à l'endroit (vue opposée)



% Puis les points cliqués

x_points_pix=-POS_objets_mesange.POS(:,2); %direction horizontale du capteur
y__points_pix=-POS_objets_mesange.POS(:,3); % direction verticale du capteur




x_0=-3840/2*ones(size(x_points_pix));
y_0=-2160/2*ones(size(x_points_pix));


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y_points = H/tanalpha_0 - H *(tanalpha_0 *(y__points_pix-y_0)+f)./(f*tanalpha_0 - (y__points_pix-y_0));

tanbeta_points=(y__points_pix-y_0)/f;
tanalpha_y_points=(tanalpha_0-tanbeta_points)./(1+ tanalpha_0*tanbeta_points);

X_points= (x_points_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y_points.^2);

% Pas de remise à l'endroit, vue opposée



% Remise en 0,0 de la caisse noire. Sur l'image, sa position est en X_points(1),Y_points(1) 
% On fait
X=X-X_points(1);
Y=Y-Y_points(1);
  
% Recalage des points à l'origine de la caisse noire
X_points=X_points-X_points(1);
Y_points=Y_points-Y_points(1);


% Rotation à corriger

% tel1  
X1=X_points(17); 
Y1=Y_points(17);
% tel2  
X2=X_points(14); 
Y2=Y_points(14);


psi = atan ((Y2-Y1)/(X2-X1));


X_prime=X;
Y_prime=Y;

X = cos(psi)*X_prime + sin(psi) *Y_prime;
Y = -sin(psi)*X_prime + cos(psi) *Y_prime;

clear X_prime Y_prime

X_points_prime=X_points;
Y_points_prime=Y_points;


X_points = cos(psi)*X_points_prime + sin(psi) *Y_points_prime;
Y_points = -sin(psi)*X_points_prime + cos(psi) *Y_points_prime;

clear X_points_prime Y__points_prime

% Figure de l'image redressée


figure(1)
surf(X,Y,image_mesange(:,:,1))
Z1=max(image_mesange(:,:,1),[],'all');

colormap(gray)
shading interp
axis equal
view(2)
title('mesange')
hold on


% On les ajoute sur le plot
%figure(1)
plot3(X_points,Y_points,1.1*double(Z1)*ones(size(X_points)),'rd','MarkerSize',18,'MarkerFaceColor','r');


figure(2)
xlim([min(min(X)) max(max(X))]);
ylim([min(min(Y)) max(max(Y))]);

%% Export to save

if ~isfolder(num2str(f))

mkdir(num2str(f))


saveas(gcf,strcat(num2str(f),'/bernache.png'))

figure(1)

saveas(gcf,strcat(num2str(f),'/mesange.png'))

end

end