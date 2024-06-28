close all
clear all


%% Redressement de l'image et passage en coordonnées physique


k=rdir('Param*');

load(k{1});

a=imread('Debut_PIV_POS_DJI_20240211153234_0273_D.tif');

x_pix=repmat([size(a,2):-1:1],size(a,1),1);
y_pix=repmat([size(a,1):-1:1]',1,size(a,2));


%x_pix=-POS(:,2); %direction horizontale du capteur
%y_pix=-POS(:,3); % direction verticale du capteur




f=2830;%f=2830; %Focale apparente en pixels

H=param.H; % Hauteur en mètres
alpha_0=param.alpha_0*pi/180; %angle passé en radians
tanalpha_0=tan(alpha_0);
% Les distances sont à prendre par rapport au centre du capteur

x_0=3840/2; % correction à faire, demi-pixel !!
y_0=2160/2;


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));

tanbeta=(y_pix-y_0)/f;

tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);

X= (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Remise à l'endroit

X=-X;
Y=-Y;

surf(X,Y,a(:,:,1))
Z1=max(a(:,:,1),[],'all');

colormap(gray)
shading interp
axis equal
view(2)
hold on

run('Redress_calib_inclin')

plot3(X,Y,double(Z1)*ones(size(X)),'rd','MarkerSize',18,'MarkerFaceColor','r');
