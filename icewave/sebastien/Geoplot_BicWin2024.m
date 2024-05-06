%% This script can be used to plot specific point son a geographical map 

% Coordinates of the relevant geographic points
% Ister, HaHa!, Saguenay , Quebec
lat = [48.452680,48.350493,48.253595,46.799582];
longi = [-68.510709,-68.809349,-70.134042,-71.221435];

txt = {'Rimouski' ,'Ha!Ha! Bay', 'Saguenay River'};

figure(1)
geoplot(lat,longi,'ob','MarkerFaceColor','b','MarkerSize',12)
geolimits manual
hold on
geobasemap('darkwater')
% text(lat,longi,txt)
gx = gca;
set(gca,'FontSize',14)
set_Papermode(gcf)


%% Extraction picture from movie
% Movie after 0211 experiment
base = 'E:/Rimouski_2024/Data/2024/0211/Drones/mesange/doc_001/';
filename = [base 'DJI_20240211213754_0224_D.MP4'];

v = VideoReader(filename);
%% Movie during 0221, multi-instrument
base = 'W:/Banquise/Présentations/Photographies/Couder_2024/0221_HaHa_deploiement/';
filename = [base 'DJI_20240221184539_0047_D.MP4'];

v = VideoReader(filename);


%% 

img = read(v,4200);

figure(2)
imshow(img)
set_Papermode(gcf)

%% Movie during 0306, multi-instrument

base = 'W:/Banquise/Présentations/Photographies/Couder_2024/0306/1-doc-vueensemble/';
filename = [base 'DJI_20240307114210_0031_D.MP4'] ;

v = VideoReader(filename);

%%
img = read(v,1400);

figure(3)
imshow(img)
set_Papermode(gcf)
