clear all
close all 

%% Loading data 

base = 'E:/Stage MIZ/PIVlab_drone/matData/';
filename = [base 'FFT_cut_50_addpad_3.mat'];

%%
disp('Loading Data..');
load(filename);
disp('Data loaded');
%% Creates a folder where to save the generated plots
fig_folder = 'E:/Stage MIZ/PIVlab_drone/Figures_FFT_cut/';
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Scaling 

fe = 30; % Frame rate in Hz
W = 64; % size of the window for DIC algorithm
font_size = 13;
fx_pix = 29.6717; % scale in pixels / meter
fx = fx_pix*2/W;
% fx = 0.8857; % for W = 64; 
%fx = 0.8857*2; % spatial scale in boxes/meter
Dt = 4; % step between two frames that were compared during the PIV algorithm 
scale_V = (fe/Dt) / fx_pix; % scale of the velocity in m/s

%% Plot the FFT spectrum 

save_boolean = 0;
[nx,ny,padding_length,n_bin] = size(TF_cut_t);

TF_spectrum = squeeze(mean(mean(mean(abs(TF_cut_t),4),2),1)); % average over all cuts and over space
N = padding_length;

TF_smooth = TF_spectrum(1:N/2+1);
TF_smooth(2:end-1) = 2*TF_smooth(2:end-1); % multiply by 2 for peaks that are both in positive an negative frequencies

f = fe*(0:(N/2))/N;

% plot FFT spectrum, fe needs to be correctly chosen 
fig_FFT_scaled = figure(6);
loglog(f,TF_smooth .* scale_V)
xlabel('$f$ (Hz)','Interpreter','latex')
ylabel('$\overline {\langle V_x \rangle _{x,y}} (f) \: \rm (m.s^{-1})$','Interpreter','latex')
ax = gca;
ax.FontSize = font_size;
grid on
axis([0.001 50 1/1e6 0.01])

% set correctly the image position for a pdf format 
set(fig_FFT_scaled,'Units','Inches');
pos = get(fig_FFT_scaled,'Position');
set(fig_FFT_scaled,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_FFT_scaled,'filename','-dpdf','-r0')

%fig_file = [fig_folder prefixe suffixe '_FFT'];
fig_file_2 = [fig_folder prefixe suffixe '_FFT_smooth_scaled'];
% saving graphs
if save_boolean 
%     saveas(fig_FFT, fig_file, 'fig');
%     saveas(fig_FFT, fig_file, 'pdf');
    saveas(fig_FFT_scaled, fig_file_2, 'fig');
    saveas(fig_FFT_scaled, fig_file_2, 'pdf');
end 


%% Plot FFT spectrum for different positions on the image 

