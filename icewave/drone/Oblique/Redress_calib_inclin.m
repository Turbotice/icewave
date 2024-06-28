%close all
%clear all


%% Redressement des images


k=rdir('Param*');
g=rdir('POS*.mat');

for i=1:length(g)
    
    load(k{i});
    load(g{i});

x_pix=-POS(:,2); %direction horizontale du capteur
y_pix=-POS(:,3); % direction verticale du capteur
f=2830;%f=2830; %Focale apparente en pixels

H=param.H; % Hauteur en mètres
alpha_0=param.alpha_0*pi/180; %angle passé en radians
tanalpha_0=tan(alpha_0);
% Les distances sont à prendre par rapport au centre du capteur

x_0=-3840/2*ones(size(x_pix));
y_0=-2160/2*ones(size(x_pix));


% On calcule les distances dans le plan à partir des coordonnées capteurs.

Y = H/tanalpha_0 - H *(tanalpha_0 *(y_pix-y_0)+f)./(f*tanalpha_0 - (y_pix-y_0));

tanbeta=(y_pix-y_0)/f;

tanalpha_y=(tanalpha_0-tanbeta)./(1+ tanalpha_0*tanbeta);

X= (x_pix-x_0)/f*H.*sqrt(1+1./tanalpha_y.^2);

% Remise à l'endroit

X=-X;
Y=-Y;

plot(X,Y,'rd','MarkerSize',18,'MarkerFaceColor','r');

% axis square
% axis equal
% hold on


end

