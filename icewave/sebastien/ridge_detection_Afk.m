%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing
date = '0223';
drone_ID = 'mesange';
exp_ID = '35-waves_014';
base = ['F:/Rimouski_2024/Data/' date '/' drone_ID '/matData/' exp_ID '/'];

filename = 'PIV_processed_i00_Dt4_b1_W32_xROI650_width3190_yROI1_height2159_scaled.mat';
matname = [base filename];

path2functions = 'C:/Users/sebas/git/icewave/drone/Drone_banquise_analysis/'; 
addpath(path2functions)

base_fig = base;
fig_folder = [base_fig 'Plots/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end
%% Load Data for A(f,k) plot
filename = ['Data_Afk_2024_' date '_' drone_ID '_' exp_ID '.mat'];

disp('Loading data..')
load([fig_folder filename])
disp('Data loaded')

%% Plot 

fig_Afk = figure; 
pcolor(k,omega,E)
shading interp
xlabel('$k \: \rm (m^{-1})$')
ylabel('$\omega \: \rm (rad.s^{-1})$')
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'ColorScale','log')
% colormap(slanCM('gnuplot2'))
cbar = colorbar();
cbar.Label.String = '$|\hat{V}_x|(k,\omega)$';
cbar.Label.Interpreter = 'latex';
axis([0.15 6, 2*pi*0.1 2*pi*2.2])
set_Papermode(gcf)
set(gca,'FontSize',13)
clim([8e-5 3e-2])

%% Try tfridge

fridge = tfridge(E,f);
fig_Afk = figure; 
pcolor(k,omega,E)
shading interp

hold on 
plot(k,2*pi*fridge,'r')