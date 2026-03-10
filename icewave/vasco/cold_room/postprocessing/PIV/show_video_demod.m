date = '20240604';
f_exc = 80;
freq_acq = 79.2;

%f_exc = 500;
%freq_acq = 124.69;

%list_f_exc = 

W = 64;
Dt = 16;

optional_intermediate_dir = '';
%optional_intermediate_dir_2 = 'sans_ondes';
%optional_sufix = '_0.3V_0.2A_cam_ext';
optional_sufix = '';


Data_demod = load(['X:/Banquise/Vasco/Frigo_pmmh/' date '/' optional_intermediate_dir '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/matData' optional_sufix '/video_demod_W' num2str(W) '_Dt' num2str(Dt) '/figdata_complex.mat']);
%Data_demod_2 = load(['X:/Banquise/Vasco/Frigo_pmmh/' date '/' optional_intermediate_dir_2 '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/matData' optional_sufix '/video_demod_W' num2str(W) '_Dt' num2str(Dt) '/figdata_complex.mat']);


%Data_demod = load(['X:/Banquise/Vasco/Frigo_pmmh/' date '/dataset2/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/matData/video_demod_W' num2str(W) '_Dt' num2str(Dt) '/figdata_complex.mat']);
disp('data loaded!')

complex_field = Data_demod.data;
%complex_field_2 = Data_demod_2.data;

%complex_field = complex_field_1 - conj(complex_field_2);

%disp(mean(abs(complex_field_1),'all'));
%disp(mean(abs(complex_field_2),'all'));

%% video imagsc
figure()
%imagesc(real(complex_field)./abs(complex_field));
nb_periods = 10;
for i=0:47+48*(nb_periods-1)%nb_periods*48-1
    %img = imagesc(real(complex_field.*exp(sqrt(-1)*2*pi*i/48).*ones(size(complex_field)))./max(abs(complex_field),[],'all'));
    %img = real(complex_field*exp(2*pi*sqrt(-1)*i/48))./max(abs(complex_field),[],'all');
    %img = img(:,5:25);
    img = real(complex_field*exp(-2*pi*sqrt(-1)*i/48))./abs(complex_field);
    imagesc(transpose(img))
    title(num2str(i))
    colorbar()
    pause(0.1)
end
%% on regarde sur un profil
index_profile_line = 5;


figure()
m = nb_periods*48;
%n = size(img,1);
%matrix = zeros(m,n);
for i=0:m-1
    img = real(complex_field*exp(2*pi*sqrt(-1)*i/48))./max(abs(complex_field),[],'all');
    %img = real(complex_field*exp(2*pi*sqrt(-1)*i/48))./abs(complex_field);
    array_to_plot = img(:,index_profile_line);
    %array_to_plot = img(30,5:25);
    plot(array_to_plot);
    %plot(transpose(img(:,25)));
    %plot(transpose(img(:,30)));
    ylim([-1 1])
    pause(0.1)
    %matrix(:,i+1) = array_to_plot;
end
%%

nb_periods = 3;
m = nb_periods*48; % Définir la valeur de m
n = size(squeeze(array_to_plot),1); % Définir la valeur de n
matrix = zeros(m,n);
for i = 1:m
 %   img = real(complex_field*exp(2*pi*sqrt(-1)*i/48))./abs(complex_field);
    img = real(complex_field*exp(2*pi*sqrt(-1)*i/48))./max(abs(complex_field),[],'all');
    p = img(:,index_profile_line);
    matrix(i,:) = p;
end
figure()
imagesc(matrix)
nb_fft_points = 2048;
%figure()
%matrix_fft2_shifted = fftshift(fft2(matrix,nb_fft_points,nb_fft_points));
%imagesc(abs(matrix_fft2_shifted))
[M,I] = max(matrix,[],1,"linear");
figure()
plot(I)

%% calcul de la vitesse de phase grace à la pente :

W = 32;
dcm = 5;
dpx = 554;
L_profil_pivboxunit = size(matrix,2);
disp(L_profil_pivboxunit);
L_profil_meters = L_profil_pivboxunit * (W/2) * (dcm*1e-2)/dpx;
disp(L_profil_meters);

% le 14/05 épaisseur environ 3mm
% le 15/05 épaisseur presque 1 cm 
%à modifier en fonction des donées :
%f_exc = 200;

t0_px = 80; % UTILISER LA FONCTION GETPOINTS !!!
x0_px = 30;
tf_px = 89;
xf_px = 77;

t0_sec = t0_px*(1/f_exc)/48; % 48 est le nombre de points par periode qu'on choisit pour le pas temporel
x0_m = x0_px * W/2 * (dcm*1e-2)/dpx;
tf_sec = tf_px*(1/f_exc)/48; % 48 est le nombre de points par periode qu'on choisit pour le pas temporel
xf_m = xf_px * W/2 * (dcm*1e-2)/dpx;

vphase = (xf_m-x0_m)/(tf_sec-t0_sec);
disp("vitesse de phase correspondante : "+num2str(vphase)+" m/s")
k = 2*pi*f_exc/vphase;
lambda = 2*pi/k;
disp("longueur d'onde sur le profil pour fexc=200Hz : "+num2str(lambda)+" m")

%à modifier en fonction des donées :
%f_exc = 100;

t0_px = 40;
x0_px = 76;
tf_px = 44;
xf_px = 13;

t0_sec = t0_px*(1/f_exc)/48; % 48 est le nombre de points par periode qu'on choisit pour le pas temporel
x0_m = x0_px * W/2 * (dcm*1e-2)/dpx;
tf_sec = tf_px*(1/f_exc)/48; % 48 est le nombre de points par periode qu'on choisit pour le pas temporel
xf_m = xf_px * W/2 * (dcm*1e-2)/dpx;

vphase = (xf_m-x0_m)/(tf_sec-t0_sec);
disp(vphase)
k = 2*pi*f_exc/vphase;
lambda = 2*pi/k;
disp("longueur d'onde sur le profil pour fexc=100Hz : "+num2str(lambda)+" m")



disp("DONNEES DU 15 mai : ")
%f_exc = 100;
disp('f_exc = '+num2str(f_exc))
t0_px = 40;
x0_px = 76;
tf_px = 44;
xf_px = 13;

t0_sec = t0_px*(1/f_exc)/48; % 48 est le nombre de points par periode qu'on choisit pour le pas temporel
x0_m = x0_px * W/2 * (dcm*1e-2)/dpx;
tf_sec = tf_px*(1/f_exc)/48; % 48 est le nombre de points par periode qu'on choisit pour le pas temporel
xf_m = xf_px * W/2 * (dcm*1e-2)/dpx;

vphase = (xf_m-x0_m)/(tf_sec-t0_sec);
disp("vitesse de phase correspondante : "+num2str(vphase)+" m/s")
k = 2*pi*f_exc/vphase;
lambda = 2*pi/k;
disp("longueur d'onde sur le profil pour fexc=100Hz : "+num2str(lambda)+" m")



%[value,idx_max] = max(abs(matrix_fft2_shifted),[],'all','linear');
%[row,col] = ind2sub(size(abs(matrix_fft2_shifted)),idx_max);

%k_fft_units = -(col - 1025); %* 2*pi; %attention on a multiplié par 2*pi



%k_profil = (2048/L_profil_meters) * k_fft_units;
%lambda_profil = 2*pi/k_profil;
%disp("longueur d'onde mesurée sur le profil : "+num2str(lambda_profil)+" m");

%%
