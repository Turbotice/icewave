
clear all
close all

%% Load data
base = 'W:/SagWin2024/Data/0226/Drones/mesange/matData/frac_001/Plots/';
filename = 'dispersion_relation_data_Vx_fmin0p1_fmax0p9_add_pow2.mat';

fullname = [base filename];
fig_folder = base;


%%
disp('Loading data...')
load(fullname)
disp('Data loaded')


%% Plot dispersion relation using 2D-FFT

% new_f = freq(1:end);
omega = 2*pi*freq;
wave_vector = k;

% first plot before masking
figure,
loglog(wave_vector,omega,'o');
grid on 
xlabel('$k$ $(m^{-1})$','Interpreter','latex');
ylabel('$\omega$ $(s^{-1})$','Interpreter','latex');
%% Computing filtered datas
% masking some datas
mask = (omega > 1.08) & (omega < 4.4);
fitted_omega = omega(mask);
fitted_k = wave_vector(mask);

figure,
loglog(fitted_k,fitted_omega,'o');
grid on 
xlabel('$k$ $(m^{-1})$','Interpreter','latex');
ylabel('$\omega$ $(s^{-1})$','Interpreter','latex');
%%
opposite2 = (fitted_omega < 0.95) & (fitted_k > 0.03) ;
mask2 = ~opposite2;
%mask2 = (fitted_omega > 0.59) ;

fitted_omega = fitted_omega(mask2);
fitted_k = fitted_k(mask2);

figure,
loglog(fitted_k,fitted_omega,'o');
grid on 
xlabel('$k$ $(m^{-1})$','Interpreter','latex');
ylabel('$\omega$ $(s^{-1})$','Interpreter','latex');

%%

opposite3 = (fitted_k > 0.5) & (fitted_omega < 1.2);
mask3 = ~opposite3;
fitted_omega = fitted_omega(mask3);
fitted_k = fitted_k(mask3);

figure,
loglog(fitted_k,fitted_omega,'o');
grid on 
xlabel('$k$ $(m^{-1})$','Interpreter','latex');
ylabel('$\omega$ $(s^{-1})$','Interpreter','latex');

%%
% mask4 = fitted_omega > 0.4;
% fitted_omega = fitted_omega(mask4);
% fitted_k = fitted_k(mask4);
%Saving relevant data
save_boolean = 1;

if save_boolean
    save_file_name = [fig_folder 'Disp_relation_plot.mat'];
    clear save
    disp('Loading plot data');
    save(save_file_name,'fitted_omega','fitted_k','x_bound','selected_freq','black_mask'); 
    disp('Loading plot data, DONE.');
end
%% Final Plot !
% fig_folder = '//192.168.1.70/Share/Data/0211/Drones/Fulmar/matData/SWO_FUL_20240211T203233UTC/Plots/continuous_ice/';
% data_plot = [fig_folder 'dispersion_relation_data_continuous_fmin0p15_fmax0p6_add_pow2.mat'];
load_bool = 0;
if load_bool
    load(data_plot)
    disp('Data to plot loaded');
end 

% fitted_k = k;
% fitted_omega = 2*pi*freq;

save_boolean = 1;

%Physical_parameters
g = 9.81; % intensity of gravity
h = 1.00; % water depth in meter
D = 3.2*10^7; % Flexion modulus
rho = 1000; % water density

% minimal values of the plot on x and y axis
xlim_min = 0.06;
xlim_max = 5;
ylim_min = 0.5;
ylim_max = 10;

% Computation of the different models 
k_list = linspace(xlim_min,xlim_max,100); % array of k for plots
deep_water = sqrt(g*k_list);
shallow_water = sqrt(g*h*k_list.^2);
flexural = sqrt(D*(k_list.^5)/rho);
gravito_waves = sqrt(g*k_list .* tanh(k_list*h));

% Power law fitting 
% data to plot
% mask_powerlaw = (fitted_omega < 2.0) & (fitted_k > 0.53) ;
% k_powerlaw = fitted_k(mask_powerlaw);
% omega_powerlaw = fitted_omega(mask_powerlaw);
% l1 = fminsearch(@(s)powerfit(k_powerlaw,omega_powerlaw,s),[1,1]);
% l1 = fminsearch(@(s)powerfit(fitted_k,fitted_omega,s),[1,1]);
% yth = powerfun(k_list,l1); % fitted exponential function

% Plot
relation_disp_scaled = figure(34);
%loglog(wave_vector,new_omega,'o');
loglog(fitted_k,fitted_omega,'o');
hold on 
loglog(k_list,deep_water,'k-');
hold on 
loglog(k_list,flexural,'b-');
% hold on 
% loglog(k_list,yth,'r-');
% hold on 
% loglog(k_list,gravito_waves,'k--');
% hold on 
% loglog(k_powerlaw,omega_powerlaw,'o','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor','black');
% title('Dispersion Relation');
xlabel('Wave number $k$ $(m^{-1})$','Interpreter','latex');
ylabel('Pulsation $\omega$ $(s^{-1})$','Interpreter','latex');
axis([xlim_min xlim_max ylim_min ylim_max]);
shallow_txt = ['Shallow, $h = ' sprintf('%0.1f',h) '\: \rm m$'];
% power_law_txt = ['$\omega = ' sprintf('%0.2f',l1(2)) 'k^{' sprintf('%0.2f',l1(1)) '}$'];
gravito_txt = ['$\omega = \sqrt{(g  k  tanh(k h))}$'];
elastic_txt = '$\omega = \sqrt{Dk^5/\rho}$';
legend('Data','$\omega = \sqrt{gk}$', elastic_txt, 'Interpreter','latex','Location','southeast');
grid on
set_Papermode(gcf)
set(findall(gcf,'-property','FontSize'),'FontSize',13)

% Saving the plot
%fig_folder = 'C:/Users/sebas/Stage_MIZ/Traitement_PIV_Rimouski_20230310/';
file_relation_disp_scaled = [fig_folder 'Relation_disp_2024_02_26_frac_001'];
if save_boolean 
    saveas(relation_disp_scaled,file_relation_disp_scaled,'pdf');
    saveas(relation_disp_scaled,file_relation_disp_scaled,'fig');
end



%% FUNCTION SECTION

function yth = powerfun(x,l)
    yth = l(2)*x.^l(1);
end 

function d = powerfit(x,y,l)
    yth = powerfun(x,l);
    d = sum((yth-y).^2);
end 