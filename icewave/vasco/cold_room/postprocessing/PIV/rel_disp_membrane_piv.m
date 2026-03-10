% Ce code prend en entrée un "film" de champs de vagues H.mat, générées à
% la fréquence f_exc et échantillonné à la fréquence f_acq.

close all

%% load variables
%base = 'X:/Banquise/Vasco/Frigo_pmmh/29032024/250Hz/matData/';
%base = 'X:/Banquise/Vasco/Frigo_pmmh/29032024/200Hz/matData/';
%base = 'X:/Banquise/Vasco/Frigo_pmmh/02042024/tests_avec_membrane/50Hz_a/matData/';
%date = '20240514';
%f_exc = 200;
%freq_acq = 99.50249;
%folder = ['I:/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' date '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/'];
%folder = ['X:/Banquise/Vasco/Frigo_pmmh/' date '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/'];
%base = [folder 'matData/'];
%base = 'X:/Banquise/Vasco/Frigo_pmmh/27032024/images/matData/70Hz/';
filename = 'X:/Banquise/Vasco/Frigo_pmmh/20240522/membrane/pulse/matData/PIV_processed_i00_N0_Dt1_b1_W32_full_total_processed.mat';

disp('Loading data..')
load(filename)
disp('Data loaded')
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

%% plots
%demod_fig = figure;
%imagesc(real(output)./abs(output))
%colorbar()
%caxis([0 15000])

%figname = [fig_folder 'demodulated_field'];

%saveas(demod_fig,figname,'fig')

N_frames=48;
%t=[0:N_frames-1]/freq_acq;

% cd "fig_folder"
% on sauve le tableau de données complexes associé (besoin d'en sauver un
% seul car démodulé):
%figdataname = sprintf("figdata_complex_Vx.mat");
%data = output2;
%save(fig_folder+figdataname,'data','-v7.3');
figure;

for i=30:200
    %A=repmat(phase,size(output,1),size(output,2));

    %imagesc(real(output.*A)./abs(output));
    %imagesc(m.x,flip(m.y),transpose(real(output2(xmin:xmax,ymin:ymax).*phase)./abs(output2(xmin:xmax,ymin:ymax))));
    imagesc(m.x,flip(m.y),transpose(H(:,10:25,i)));
    %imagesc(m.x,flip(m.y),transpose(real(output2.*phase)))
    caxis([-10 10])
    colorbar()
    %imwrite(real(h.*A),turbo,strcat('W:/Banquise/Samuel/champs_demodule_20230512/acq_50/','img',num2str(i),'.png'))
    daspect([1 1 1])
    %figname = sprintf("figure_%d.tiff",i);%[fig_folder 'demodulated_field_' string(i)];
    %saveas(gcf,fig_folder + figname)
    %[xi,yi] = getpts();
    %x(i) = xi;
    %y(i) = yi;
    title(num2str(i))
     pause(0.01)
    i
    %T(i) = getframe(gcf);


%     close(gcf);
end 
%%
W = 32;
%dcm = 7;
%dpx = 1192;
dcm = 11;
dpx = 1215;

matrix_spatio_temp = squeeze(H(:,21,30:900));
dx_pix = W/2;
dx_meters = ((dcm*1e-2)/dpx)*dx_pix;
fx = 1/(dx_meters);%en boites par metre
freq_acq = 169;
Nk = size(matrix_spatio_temp,1);
Nf = size(matrix_spatio_temp,2); 
freq_values = (-1024/2:1024/2-1)*freq_acq*1/1024;
k_values = (-1024/2:1024/2-1)*fx*(2*pi)*1/1024;
figure;
imagesc(k_values,freq_values,abs(fftshift(fft2(matrix_spatio_temp,1024,1024)))');
hold on
rho = 965;
e = 520e-6;
nu = 0.5;
g = 9.81;
h = 2.9e-2; % profondeur eau

E = 1.6e6;
D = (E*(e^3))/(12*(1-nu^2));
plot(k_values,(1/(2*pi))*sqrt(g*k_values+((D/rho)*k_values.^5).*tanh(h*k_values)),Color=[0 1 0]);

E = 2.45e6;
D = (E*(e^3))/(12*(1-nu^2));
plot(k_values,(1/(2*pi))*sqrt(g*k_values+((D/rho)*k_values.^5).*tanh(h*k_values)),Color=[1 0 0]);
legend('theory for E = 1.6 Mpa','theory for E = 2.45 Mpa')
xlabel('k (m^{-1})')
ylabel('f (Hz)')
set(gca,'YDir','normal')
title('Dispersion relation of hydroelastic waves')
%set(gca,'YScale','log')
%set(gca,'XScale','log')
%colorbar;
caxis([0 0.8e4]);
axis([0 500 0 83]);
%% relation de dispersion theorique
k_th = linspace(0,250,1000);
disp(D);
omega_th = sqrt((g*k_th + (D/rho)*(k_th.^5)).*tanh(k_th*h));
figure;
plot(k_th,(1/(2*pi))*omega_th);
hold on
plot(k_th,(1/(2*pi))*sqrt(10*k_th))
