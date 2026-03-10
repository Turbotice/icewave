% Ce code prend en entrée un "film" de champs de vagues H.mat, générées à
% la fréquence f_exc et échantillonné à la fréquence f_acq.

close all

%% load variables
%base = 'X:/Banquise/Vasco/Frigo_pmmh/29032024/250Hz/matData/';
%base = 'X:/Banquise/Vasco/Frigo_pmmh/29032024/200Hz/matData/';
%base = 'X:/Banquise/Vasco/Frigo_pmmh/02042024/tests_avec_membrane/50Hz_a/matData/';
date = '20240514';
f_exc = 200;
freq_acq = 99.50249;
%folder = ['I:/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' date '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/'];
folder = ['X:/Banquise/Vasco/Frigo_pmmh/' date '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/'];
base = [folder 'matData/'];
%base = 'X:/Banquise/Vasco/Frigo_pmmh/27032024/images/matData/70Hz/';
filename = [base 'PIV_processed_i00_N0_Dt8_b1_W32_full_total_processed.mat'];

disp('Loading data..')
load(filename)
disp('Data loaded')

%%
%fig_folder = 'X:/Banquise/Vasco/Frigo_pmmh/29032024/100Hz/video_demod/';
fig_folder = [folder '/matData/video_demod/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end
disp('video folder exists')
%% Création de la projection radiale
% definition des coordonnées de la source:
%xs = 0 ;
%ys = 0 ;
%xs=0;
%ys=0;
%xs = 50 ;
%ys = -10 ;
%[Y,X] = meshgrid(m.y,m.x);
%Ux = (1./sqrt((X-xs).^2 + (Y-ys).^2)) .* (X - xs);
%Uy = (1./sqrt((X-xs).^2 + (Y-ys).^2)) .* (Y - ys);
%Ux3d = repmat(Ux,[1 1 size(m.Vx,3)]);
%Uy3d = repmat(Uy,[1 1 size(m.Vy,3)]);
% Le signal d'intéret est le suivant: (projection radiales des vitesses)
%H = Ux3d .* m.Vx + Uy3d .* m.Vy;
H = m.Vy;
%H = m.Vx;

%mean_along_dim_1_2 = mean(mean(H,1),2);
%H_means = repmat(mean_along_dim_1_2, size(H, 1), size(H, 2), 1);
%H_corrected = H - H_means;
%H = H_corrected;
disp('H array created')

%%
%load('PIV_processed_i00_N0_Dt1_b1_W16_full_total_processed');
%f_exc=100;
%freq_acq=99.0099;
%f_exc = 50;
%freq_acq = 49.5049;
%f_exc = 70;
%freq_acq = 169;
%f_exc = 100;
%freq_acq = 99.0099;
%f_exc = 150;
%freq_acq = 148.5222;
% Le champ H contient 3 dimensions : H(x,y,t)
%H=m.Vy;%sqrt(m.Vx.^2+m.Vy.^2);
%H=m.Vx;%sqrt(m.Vx.^2+m.Vy.^2);


% Etape 1 : initialisations

Y=zeros(size(H,1),size(H,2));
Nframes=size(H,3);

% Vecteur temps



t=[0:Nframes-1]/freq_acq;

phase=exp(2*sqrt(-1)*pi*f_exc*t);

% Construction de la matrice de demodulation
clear A;
A=repmat(phase,size(H,2),1);
B(1,:,:)=A;
demod=repmat(B,size(H,1),1,1);
clear A;
clear B;

% On fait ensuite le produit entre H et demod (terme à terme), et la somme
% dans la direction temporelle pour sortir le coef de fourier
disp('demodulation en cours');
output=sum(H.*demod,3);
%output2 = output-mean(output,'all');
disp('demodulation teminée')

%% plots
%demod_fig = figure;
%imagesc(real(output)./abs(output))
%colorbar()
%caxis([0 15000])

%figname = [fig_folder 'demodulated_field'];

%saveas(demod_fig,figname,'fig')

N_frames=48;
t=[0:N_frames-1]/freq_acq;

nb_phases = 10;
% cd "fig_folder"
% on sauve le tableau de données complexes associé (besoin d'en sauver un
% seul car démodulé):
%figdataname = sprintf("figdata_complex_Vx.mat");
%data = output2;
%save(fig_folder+figdataname,'data','-v7.3');
    figure;

    %x = zeros(1,N_frames);
    %y = zeros(1,N_frames);

for i=1:nb_phases*N_frames
    phase=exp(-2*sqrt(-1)*pi*i/N_frames);
    %A=repmat(phase,size(output,1),size(output,2));

    %imagesc(real(output.*A)./abs(output));
    %imagesc(m.x,flip(m.y),transpose(real(output2(xmin:xmax,ymin:ymax).*phase)./abs(output2(xmin:xmax,ymin:ymax))));
    imagesc(m.x,flip(m.y),transpose(real(output.*phase)./abs(output)));
    %imagesc(m.x,flip(m.y),transpose(real(output2.*phase)))
    %caxis([-1000 1000])
    caxis([-1 1])
    colorbar()
    %imwrite(real(h.*A),turbo,strcat('W:\Banquise\Samuel\champs_demodule_20230512\acq_50\','img',num2str(i),'.png'))
    daspect([1 1 1])
    figname = sprintf("figure_%d.tiff",i);%[fig_folder 'demodulated_field_' string(i)];
    %saveas(gcf,fig_folder + figname)
    %[xi,yi] = getpts();
    %x(i) = xi;
    %y(i) = yi;

     pause(0.01)
    i
    %T(i) = getframe(gcf);


%     close(gcf);
end 

%%
%vi = VideoWriter('Video_test','avi');
%vi.FrameRate = 10;
%open(vi)
%writevideo(vi,T)
%close(vi)

%% fft2
% first choose in which region to perform the fft
xmin=40;
xmax=130;
ymin=1;
ymax=100;

w = 16;
dx_photo = 2e-4;
dx_piv = dx_photo*w/2;
% remarque : nextpow2 utile pour padding
n_pad_x=1000;
n_pad_y=1000;
figure;
FT2 = fft2(output2(xmin:xmax,ymin:ymax),size(H,1)+n_pad_x,size(H,2)+n_pad_y);
kx=linspace(-pi/dx_piv,pi/dx_piv, size(H,1)+n_pad_x);
ky=linspace(-pi/dx_piv,pi/dx_piv, size(H,2)+n_pad_y);

imagesc(kx,flip(ky),transpose(abs(fftshift(FT2))));
colorbar()
