%close all
%clear all


%% Redressement de l'image et passage en coordonnées physique


k=rdir('Param*');

load(k{1});

a=imread('Debut_PIV_DJI_20240211213233_0222_D.tif');

x_pix=repmat([size(a,2):-1:1],size(a,1),1);
y_pix=repmat([size(a,1):-1:1]',1,size(a,2));


%x_pix=-POS(:,2); %direction horizontale du capteur
%y_pix=-POS(:,3); % direction verticale du capteur




f=2830;%f=2830; %Focale apparente en pixels

H=param.H; % Hauteur en mètres
alpha_0=param.alpha_0*pi/180; %angle passé en radians
tanalpha_0=tan(alpha_0);
% Les distances sont à prendre par rapport au centre du capteur

x_0=3840/2;
y_0=2160/2;


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));

tanbeta=(y_pix-y_0)/f;

tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);

X= (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Remise à l'endroit. On ne change pas le signe pour s'aligner avec les
% axes proposés par bernache qui vont servir de ref.
X=X;
Y=Y;

% On doit caler le 0,0 sur la caisse noire

X=X+18.9588;
Y=Y+10.4427;

% Reste à tourner les axes d'un angle donné. On prend les tels de
% l'alignement de bouées et on sort leurs coordonnées avant rotation pour
% déterminer l'angle de rotation psi
% tel1  
X1=-20.1810; 
Y1= -2.0211;
% tel2  
X2= 70.0488; 
Y2= 3.0612;



psi = atan ((Y2-Y1)/(X2-X1));

% On doit maintenant appliquer une matrice de rotation d'un angle -psi sur
% X,Y

X_prime=X;
Y_prime=Y;

X = cos(psi)*X_prime + sin(psi) *Y_prime;
Y = -sin(psi)*X_prime + cos(psi) *Y_prime;





figure(1)
surf(X,Y,a(:,:,1))
Z1=max(a(:,:,1),[],'all');

colormap(gray)
shading interp
axis equal
view(2)
hold on

run('Redress_calib_inclin')
% 
plot3(X,Y,double(Z1)*ones(size(X)),'rd','MarkerSize',18,'MarkerFaceColor','r');
xlim([-80 120])
ylim([-55 45])