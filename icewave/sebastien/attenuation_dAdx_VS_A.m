%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing
date = '0223';
drone_ID = 'bernache';
exp_ID = '07-waves_007';
base = ['K:/Share_hublot/Data/' date '/Drones/' drone_ID '/matData/' exp_ID '/'];

filename = 'PIV_processed_i00_N0_Dt5_b1_W32_xROI1013_width2827_yROI1_height2159_scaled.mat';
matname = [base filename];

path2functions = 'C:/Users/sebas/git/icewave/drone/Drone_banquise_analysis/'; 
addpath(path2functions)
path2formula_functions = 'C:/Users/sebas/git/icewave/drone/Drone_banquise_analysis/formula_fit/';
addpath(path2formula_functions)

%% Load data 
disp('Loading Data..');
load(matname);
disp('Data loaded');

%%  Folder where plots are saved and set parameters for plotting 
base_fig = base;
fig_folder = [base_fig 'Plots/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end
%% Scales

facq_x = 1/m.SCALE.fx; % scale in box / meter
fx = m.SCALE.fx;
facq_t = 1/m.SCALE.ft; % scale in frame / sec
scale_V = m.SCALE.scale_V; % factor scaling for velocity in meter / s
W = m.PIV_param.w ;

%% Get rid off quadratic noise (drone movements)
[nx,ny,nt] = size(m.Vx);
x = (1:1:nx);
y = (ny:-1:1);

% compute the mean component of the velocity field for each frame
Vxmoy = mean(mean(m.Vx,2),1);
Vymoy = mean(mean(m.Vy,2),1);

Vx = supress_quadratic_noise(m.Vx,x,y); 
Vy = supress_quadratic_noise(m.Vy,x,y);
% Vx = flip(Vx,2);
% Vy = flip(Vy,2);

% Vx = flip(Vx,1); % waves are initially coming from the right boarder of the image
% Vy = flip(Vy,1);

disp('Drone motion corrected')


% #######################################
%% Attenuation coefficient in corrected direction 
% #######################################

% Perform time FFT 
disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;
[FFT_t,~,f] = temporal_FFT(Vx(:,:,:),padding_bool,add_pow2,facq_t);


% Parameters 
selected_freq = [0.25 0.6]; % selected frequencies between which we proceed to the analysis
% get indices of frequencies closest to max and min selected frequency
[min_freq, start_freq_idx] = min(abs(f - selected_freq(1))); 
[max_freq, end_freq_idx] = min(abs(f - selected_freq(2)));
f_cropped = f(start_freq_idx : end_freq_idx); % new frequency array 
FFT_cropped = FFT_t(:,:,start_freq_idx:end_freq_idx);

x = m.x;
y = m.y;

L0 = 100; % segments size over which we look at attenuation
xmin = 1; % minimal value of x-coordinate over which we fit attenuation 

% Select range of frequencies and space of shortened fit 
f_cutoff = [0.43 0.46 0.49 0.53 0.57]; % frequency above which the attenuation fit is shortened
x_cutoff = [80 60 45 35 25]; % meters % range over which the shortened attenuation fit is performed

new_folder_fig = [fig_folder 'correct_orientation_A_VS_x/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

[S_A] = amplitude_corrected_direction(FFT_cropped,...
    f_cropped,x,y,facq_x,L0,xmin,f_cutoff,x_cutoff,new_folder_fig);

%% Create a structure and save it

S_Ax = struct('S_A',S_A,'f',f_cropped,'L0',L0,'f_cutoff',f_cutoff,'x_cutoff',x_cutoff);
 
attenuation_file =  ['Data_Ax_' date '_' drone_ID '_' exp_ID '_' num2str(selected_freq(1)) 'Hz_to_' num2str(selected_freq(2)) 'Hz'];
attenuation_file = replace(attenuation_file,'.','p');
attenuation_file = [fig_folder attenuation_file '_harmonic_1'];
save(attenuation_file,'S_Ax')

disp('DONE.')


%% Fit each profile A(x) by an exponential 

nf = length(S_Ax.S_A);
alpha = zeros(nf,1); % array of attenuation coefficients
C  = zeros(nf,1); % array of prefactor
d = zeros(nf,1); % array of distance to plot

for idx_freq = 1 : length(S_Ax.S_A)
   
    A_red = S_Ax.S_A(idx_freq).A;
    x_fit = S_Ax.S_A(idx_freq).x;

    [alpha,C,d] = exponential_coeff(A_red',x_fit);
    alpha(idx_freq) = alpha;
    C(idx_freq) = C;
    d(idx_freq) = C;

end 

%% Save attenuation structure 

S_atten = struct('alpha',alpha,'C',C,'d',d,'f',S_Ax.f,'L0',L0,'f_cutoff',f_cutoff,'x_cutoff',x_cutoff);
 
attenuation_file =  ['Data_attenuation_' date '_' drone_ID '_' exp_ID '_' num2str(selected_freq(1)) 'Hz_to_' num2str(selected_freq(2)) 'Hz'];
attenuation_file = replace(attenuation_file,'.','p');
attenuation_file = [fig_folder attenuation_file '_harmonic_1'];
save(attenuation_file,'S_atten')

disp('DONE.')

%% Plot distance d to the curve as function of frequencies
figure, 

plot(S_atten.f,S_atten.d,'o')
xlabel('$f \: \rm (Hz)$')
ylabel('$d$')


%% 

d_thresh = 0.05;

mask = (S_atten.d < d_thresh) & (abs(S_atten.alpha) > 0.007) & (S_atten.f' < 0.6) & (S_atten.f' > 0.3);
    % (S_atten.f' < 0.55) & (S_atten.f' > 0.18); % keep only points for which distance is smaller than..
freq = S_atten.f;
fitted_f = freq(mask);
fitted_alpha = abs(S_atten.alpha(mask))';

l1 = fminsearch(@(s)powerfit(fitted_f,fitted_alpha,s),[1,1]);
f_list = linspace(0.01,10,100);
yth = powerfun(f_list,l1); % fitted powerlaw function

attenuation_fig = figure;
loglog(fitted_f,fitted_alpha,'o','MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerEdgeColor','black');
hold on
plot(f_list,yth,'r--','LineWidth',1.5);
xlabel('$f \: \rm (Hz)$','Interpreter','latex');
ylabel('$\alpha \: \rm (m^{-1})$','Interpreter','latex');
grid on 
axis([1e-1 2 1e-3 1])
ax = gca;
ax.FontSize = 13;

power_law_txt = ['$\alpha(f) = ' sprintf('%0.2f',l1(2)) 'f^{' sprintf('%0.2f',l1(1)) '}$'];
% legend('Data','Fitted Data',power_law_txt,'Interpreter','latex','Location','southeast','FontSize',13)
legend('',power_law_txt,'Interpreter','latex','Location','northwest','FontSize',13)
set_Papermode(gcf);

thresh_txt = replace(num2str(1-d_thresh),'.','p');
attenuation_filename = [fig_folder 'attenuation_law_' date '_' drone_ID '_' exp_ID 'harmonic_1_confidence_' thresh_txt];
saveas(attenuation_fig,attenuation_filename,'fig');
saveas(attenuation_fig,attenuation_filename,'pdf');

%% Load amplitudes

attenuation_file =  ['Data_Ax_' date '_' drone_ID '_' exp_ID '_' num2str(selected_freq(1)) 'Hz_to_' num2str(selected_freq(2)) 'Hz'];
attenuation_file = replace(attenuation_file,'.','p');
attenuation_file = [fig_folder attenuation_file '_harmonic_1.mat'];

if isfile(attenuation_file)
    disp('Attenuation file already exists')
    load(attenuation_file);
    disp('Attenuation data loaded')
else
    disp('No attenuation file saved')
end 


%% Differentiate a profile A(x)
S_A = S_Ax.S_A;

f0 = 0.42;
[~,i0] = min(abs(S_Ax.f - f0));

% Parameters for Savitzky-Golay method
order = 3;
ratio = 3;


% Set sub_fig_folder
sub_fig_folder = [fig_folder 'dAdx/'];
if exist(sub_fig_folder,'dir') ~= 7
    mkdir(sub_fig_folder)
end

for i0 = 1: length(S_Ax.f)
    disp(['f = ' num2str(S_Ax.f(i0),'%.4f') ' Hz'])
    freq_txt = replace(num2str(S_Ax.f(i0)),'.','p'); % suffixe used to name figures

    framelen = 2*round(length(S_A(i0).A)/ratio/2) + 1;
    disp(framelen)
    sgf = sgolayfilt(S_A(i0).A,order,framelen);
    
    profile_fig = figure(6);
    plot(S_A(i0).x,S_A(i0).A,'-')
    hold on
    plot(S_A(i0).x,sgf,'.-')
    
    % differentiate
    dA = gradient(S_A(i0).A,1/facq_t);
    dsgf = gradient(sgf,1/facq_t);
    
    % fit by exponential
    [alpha,C,d] = exponential_coeff(sgf',S_A(i0).x);
    y_th = C*exp(alpha.*S_A(i0).x);
    
    hold on
    plot(S_A(i0).x,-dsgf,'.-')
    hold on 
    plot(S_A(i0).x,y_th,'--')
    
    xlabel('$x$')
    ylabel('$A(x)$')
    legend('signal','sgolay','dsgolay','exp')
    hold off

    % save figure
    figname = [sub_fig_folder 'Ax_profile_' freq_txt];
    saveas(profile_fig,figname,'fig')
    saveas(profile_fig,figname,'pdf')
    hold off

    dAdx_fig = figure(7);
    dexp = gradient(y_th,1/facq_t);
    
    plot(sgf,-dsgf,'.-','Color',[0.8500 0.3250 0.0980])
    hold on 
    plot(y_th,-dexp,'--','Color',[0.4940 0.1840 0.5560])
    xlabel('$A$')
    ylabel('$- \frac{dA}{dx}$')
    legend('sgf','exponential')


    figname = [sub_fig_folder 'dAdx_' freq_txt];
    saveas(dAdx_fig,figname,'fig')
    saveas(dAdx_fig,figname,'pdf')
    hold off
end



% ##########################
%% Load image 

path2img = ['K:/Share_hublot/Data/' date '/Drones/' drone_ID '/' exp_ID '/*exemple.tiff'];
filelist = dir(path2img);

path = [filelist(1).folder '/' filelist(1).name];
img = imread(path);

%%

img_dim = size(img);
yWorldLimits = [0,img_dim(1)-1]/m.SCALE.facq_pix;
xWorldLimits = [0,img_dim(2)-1]/m.SCALE.facq_pix;
R = imref2d(img_dim(1:2),xWorldLimits,yWorldLimits);

figure(8)
imshow(img,R)
set(gca,'YDir','normal')

%%
cropped_img = img(m.PIV_param.ROI.y : m.PIV_param.ROI.y + m.PIV_param.ROI.height,...
                  m.PIV_param.ROI.x : m.PIV_param.ROI.x + m.PIV_param.ROI.width);

img_dim = size(cropped_img);
yWorldLimits = [0,img_dim(1)-1]/m.SCALE.facq_pix;
xWorldLimits = [0,img_dim(2)-1]/m.SCALE.facq_pix;
R = imref2d(img_dim(1:2),xWorldLimits,yWorldLimits);

fig_cropped = figure(8);
imshow(cropped_img,R)
set(gca,'YDir','normal')

figname = [fig_folder 'Area_PIV_scaled_' drone_ID '_' exp_ID];
saveas(profile_fig,figname,'fig')
saveas(profile_fig,figname,'pdf')
hold off

%% Create a function to scale images 

%% Try exponential fit using two parameters ?? 








