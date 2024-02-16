
clear all
close all

%% Load datas Diagonal Scaling

% Load datas associated to waves
folder_waves = 'E:/Stage MIZ/PIVlab_drone/Figures/DJI_0402_figures_Dt2_W32_waves_disp_relation/';
waves_data = [folder_waves 'Relevant_datas_waves_2023_11_22_diagonal_scaling_addpad3.mat'];
load_bool = 1;
if load_bool
    waves = load(waves_data); % store fitted_k and fitted_omega in a structure "waves"
    disp('Data to plot loaded');
end 

% Load datas associated to fragmented ice
folder_ice = 'E:/Stage MIZ/PIVlab_drone/Figures/DJI_0402_figures_Dt4_disp_relation/';
ice_data = [folder_ice 'Relevant_datas_ice_2023_11_22_diagonal_scaling_addpad2.mat'];
load_bool = 1;
if load_bool
    ice = load(ice_data); % store fitted_k and fitted_omega in a structure "ice"
    disp('Data to plot loaded');
end 

fig_folder = 'E:/Stage MIZ/PIVlab_drone/Figures_report/';
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Load datas Scaling using theta_x

% Load datas associated to waves
folder_waves = 'E:/Stage MIZ/PIVlab_drone/Figures_report/new_scale_waves/';
waves_data = [folder_waves 'Relevant_datas_waves_2023_12_04_Elie_scaling_addpad3.mat'];
load_bool = 1;
if load_bool
    waves = load(waves_data); % store fitted_k and fitted_omega in a structure "waves"
    disp('Data to plot loaded');
end 

% Load datas associated to fragmented ice
folder_ice = 'E:/Stage MIZ/PIVlab_drone/Figures_report/new_scale_ice/';
ice_data = [folder_ice 'Relevant_datas_ice_2023_12_04_Elie_scaling_addpad3.mat'];
load_bool = 1;
if load_bool
    ice = load(ice_data); % store fitted_k and fitted_omega in a structure "ice"
    disp('Data to plot loaded');
end 

fig_folder = 'E:/Stage MIZ/PIVlab_drone/Figures_report/';
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Plot 
save_boolean = 0;

%Physical_parameters
g = 9.81; % intensity of gravity
h = 3.9; % water depth in meter
D = 3.2*10^7; % Flexion modulus
rho = 1000; % water density

% minimal values of the plot on x and y axis
xlim_min = 0.03;
xlim_max = 10;
ylim_min = 0.4;
ylim_max = 7;

% Computation of the different models 
k_list = linspace(xlim_min,xlim_max,100); % array of k for plots
deep_water = sqrt(g*k_list);
shallow_water = sqrt(g*h*k_list.^2);
flexural = sqrt(D*(k_list.^5)/rho);
gravito_waves = sqrt(g*k_list .* tanh(k_list*h));

% Power law fitting 
% data to plot
% mask_powerlaw = (fitted_omega < 2.0) & (fitted_k > 0.6) ;
% k_powerlaw = fitted_k(mask_powerlaw);
% omega_powerlaw = fitted_omega(mask_powerlaw);
% l1 = fminsearch(@(s)powerfit(k_powerlaw,omega_powerlaw,s),[1,1]);
% % l1 = fminsearch(@(s)powerfit(fitted_k,fitted_omega,s),[1,1]);
% yth = powerfun(k_list,l1); % fitted exponential function

% Plot
relation_disp_scaled = figure(1);
%loglog(wave_vector,new_omega,'o');

loglog(waves.fitted_k,waves.fitted_omega,'o','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor','black','MarkerSize',6);
hold on 
loglog(ice.fitted_k,ice.fitted_omega,'o','MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerEdgeColor','black');
hold on 
loglog(k_list,deep_water,'k--','LineWidth',1.5);
hold on 
loglog(k_list,shallow_water,'--','LineWidth',1.5,'Color','r');
% hold on 
% loglog(k_list,gravito_waves,'k--');
%title('Dispersion Relation','Interpreter','latex');
xlabel('$k$ $(\rm m^{-1})$','Interpreter','latex');
ylabel('$\omega$ $(\rm s^{-1})$','Interpreter','latex');
axis([0.04 10 0.1 10]);
deep_txt = ['$\omega = \sqrt{(g  k )}$'];
shallow_txt = ['$\omega = \sqrt{(g  h_w k^2)}$, $h_w = ' sprintf('%0.1f',h) '\: \rm m$'];
%power_law_txt = ['$\omega = ' sprintf('%0.2f',l1(2)) 'k^{' sprintf('%0.2f',l1(1)) '}$'];
gravito_txt = ['$\omega = \sqrt{(g  k \; \tanh(k h_w))}$'];
legend('Waves','Ice',deep_txt, shallow_txt, 'Interpreter','latex','Location','northwest','location','southeast');
grid on

ax = gca;
ax.FontSize = 13;
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'filename','-dpdf','-r0')

% Saving the plot

file_relation_disp_scaled = [fig_folder 'Relation_disp_scaled_complete_2023_12_04'];
if save_boolean 
    saveas(relation_disp_scaled,file_relation_disp_scaled,'pdf');
    saveas(relation_disp_scaled,file_relation_disp_scaled,'fig');
end