

%% Load data 
fps = 100; % frame rate
facq_x = 2400; % pix/meter

base = 'F:/Waves_reconstruction_wilson/9_wave_field_2024_01_02/';
directory = [base 'D4cm_h10mm_fps_100/f6.0Hz_amplitude15mm/'];

crop_matrix_file = [directory 'croppedMatrixobject.mat'];
load(crop_matrix_file)

pos_file = [directory 'locss.mat'];
load(pos_file)

freq_file = [directory 'freqdemod.mat'];
load(freq_file)

%% Choosing the position of the profile in the matrix using locs
[nx,ny,nt]=size(croppedMatrix);
first_frame=1;
last_frame=nt;
fx = facq_x;

% load croppedMatrix.mat croppedMatrix
% load locss.mat locs
prof = round(min(locs(:,1))-10); % change this line !!

%% taking profile
% [nx,ny,nt]=size(croppedMatrix);
first_frame=1;
last_frame=nt;

freq_acq = fps;
f_exc_range = f_subpix;  % f_exc value for demodulation

radious= 2; %cm

for i=first_frame:last_frame
 wPRc_v(:,1,i)=croppedMatrix(:,prof, i);
end

%% FOURIER ANALYSIS AND FILTERING
ke = 1./fx; %% echantillonage en k
m.k=(-ke/2:ke/(nx-1):ke/2); % axe des k
m.f =(-fps/2:fps/(nt-1):fps/2); % axe des temps
trv=squeeze(wPRc_v); % 1st dimension is space, 2nd is time
figure(31)
subplot(3,3,1)
imagesc(trv)
%savefig('spacetemporal.fig') 
tr1v=fftshift(fft2(trv)); % Computes the FFT both in space and time for the spaciotemporal
% tr1v=tr1v;

subplot(3,3,2)
imagesc(((abs(tr1v))))
% imagesc(KK,KK,(log(abs(tr1v))))
 ylabel('px')
 xlabel('frames')

subplot(3,3,3)
 imagesc(m.f,m.k,(abs(tr1v)))
 
 ylabel('kx(cm-1)')
 xlabel('f(Hz)')
 xlim([-6 6])
 ylim([-1 1])
 
 subplot(3,3,4) %Original
imagesc(trv)

subplot(3,3,5)%Isolated

Y_v=tr1v;

% Y_v(248.5:252.5,1:nt)=zeros;
% Y_v(1:end,498:504)=zeros;

%TO KEEP WAVES MOVING TO THE RIGHT of the cylinder
Y_v(1:round(nx/2),1:round(nt/2))=zeros;
Y_v(round(nx/2):end,(round(nt/2)):nt)=zeros;

%BAND high K
% Y_v(1:240,1:end)=zeros;
% Y_v(260:end,1:end)=zeros;

signal_disk_v1=real(ifft2((ifftshift(Y_v))));
imagesc((signal_disk_v1))

wave2=tr1v-Y_v;

subplot(3,3,6)%Result

signal_wavev1=real(ifft2(ifftshift(wave2)));

imagesc((signal_wavev1))
 
subplot(3,3,7) %Original
imagesc(trv)

subplot(3,3,8)%Isolated
Y_v=tr1v;

%TO KEEP WAVES MOVING TO THE LEFT  of the cylinder
Y_v(1:round(nx/2),round(nt/2):end)=zeros;
Y_v(round(nx/2):end,1:(round(nt/2)))=zeros;

%BAND high K
% Y_v(1:240,1:end)=zeros;
% Y_v(260:end,1:end)=zeros;

signal_disk_v2=real(ifft2((ifftshift(Y_v))));
imagesc((signal_disk_v2))

wave2=tr1v-Y_v;

subplot(3,3,9)%Result

signal_wavev2=real(ifft2(ifftshift(wave2)));

imagesc((signal_wavev2))
 %savefig('fft.fig') 

 trvresult1 = permute(signal_disk_v1, [3, 1, 2]);
 trvresult2 = permute(signal_disk_v2, [3, 1, 2]);
 
rigth=trvresult1;  %isolated signal going from the cylinder to the rigth
left=trvresult2;   %isolated signal going from the cylinder to the left

%% putting the cylinder at the center of the axis
 d = size(trv);
%  X = linspace(0,max(d),max(d));
 Y = linspace(0,min(d),min(d));
s=Y(1,1:nx)*fx;
max_val = max(s(:));
% calculate the midpoint
if mod(max_val,2) == 0 % even number of elements
    midpoint = max_val/2 + 0.5; % add 0.5 to shift midpoint between two elements
else % odd number of elements
    midpoint = (max_val/2);
end

% subtract the midpoint from each element
shifted_matrix = s - midpoint;
s1=shifted_matrix;


%% Join left and rigth signal 

trv_p = permute(trv, [3, 1, 2]);


[nx,ny,nt]=size(trv_p);
first_frame=1;
last_frame=nt;
% fps=100;
% dcm = 6 % cm rule on the surface of the liquid
% dpx = 143 % px
% fx  = dcm/dpx % ratio cm/px
% scale=dpx/(dcm*10)%15; % pixels per mm
% vertical


for i=1:last_frame%length(g)

indicesb = abs(shifted_matrix) <= radious; 
    trv_p(1 ,indicesb, i) = NaN;
end
 
for i=1:last_frame%length(g)

indicesr = shifted_matrix <= radious; 
    rigth( 1,indicesr, i) = 0;
 
indicesl = shifted_matrix >= -radious;   
    left( 1,indicesl, i) = 0;

end

trv_new=rigth+left;  

trv_newtest=trv_new;
        

s = {}
for i=1:last_frame%length(g)0 


indicesr = abs(shifted_matrix) <= (radious+0.1); 
    trv_newtest( 1,indicesr, i) = NaN;
    
    % Find the indices of NaN elements
nanIndices( 1,:, 1) = isnan(trv_new( 1,:, 1));
s(i).trv_new = trv_new(1,:,i);

NONaN( 1,:, i) = s(i).trv_new(~nanIndices);


end

%% Band pass filter to get rid off low and high wavenumbers
for i=1:nt
   
    trvBP(:,1,i) = filt_band_pass2(trv_new(:,:,i),[5 30],0); %apply band pass filter
   disp(i)
end


%% Visualization step
figure(14)
ke = 1./fx; %% echantillonage en k
m.k=(-ke/2:ke/(ny-1):ke/2); % axe des k


m.f =(-fps/2:fps/(nt-1):fps/2); % axe des temps


% trv=squeeze(wPRc_v); % 1st dimension is space, 2nd is time
% figure(31)
subplot(3,4,1)
imagesc(trv)
%savefig('spacetemporal.fig') 
tr1v=fftshift(fft2(trv)); % Computes the FFT both in space and time for the spaciotemporal
% tr1v=tr1v;

subplot(3,4,2)
imagesc(((abs(tr1v))))
% imagesc(KK,KK,(log(abs(tr1v))))
 ylabel('px')
 xlabel('frames')

subplot(3,4,3)
 imagesc(m.f,m.k,(abs(tr1v)))
 
 ylabel('kx(cm-1)')
 xlabel('f(Hz)')
 xlim([-10 10])
 ylim([-1 1])

 
subplot(3,4,4)%Resultmage
imagesc(m.f,m.k,log(abs(tr1v)))
 ylabel('kx(cm-1)')
 xlabel('f(Hz)')
 xlim([-10 10])
 ylim([-1 1])

 
 subplot(3,4,5)
trp=squeeze(trv_new);
imagesc(trp)
%savefig('spacetemporal.fig') 
tr1p=fftshift(fft2(trp)); % Computes the FFT both in space and time for the spaciotemporal
% tr1v=tr1v;

subplot(3,4,6)
imagesc(((abs(tr1p))))
% imagesc(KK,KK,(log(abs(tr1v))))
 ylabel('px')
 xlabel('frames')

subplot(3,4,7)
 imagesc(m.f,m.k,(abs(tr1p)))
 
 ylabel('kx(cm-1)')
 xlabel('f(Hz)')
 xlim([-10 10])
 ylim([-1 1])

 subplot(3,4,8)%Resultmage
imagesc(m.f,m.k,log(abs(tr1p)))
 ylabel('kx(cm-1)')
 xlabel('f(Hz)')
 xlim([-10 10])
 ylim([-1 1])
 
  subplot(3,4,9)
trb=squeeze(trvBP);
imagesc(trb)
%savefig('spacetemporal.fig') 
tr1b=fftshift(fft2(trb)); % Computes the FFT both in space and time for the spaciotemporal
% tr1v=tr1v;

subplot(3,4,10)
imagesc(((abs(tr1b))))
% imagesc(KK,KK,(log(abs(tr1v))))
 ylabel('px')
 xlabel('frames')

subplot(3,4,11)
 imagesc(m.f,m.k,(abs(tr1b)))
 
 ylabel('kx(cm-1)')
 xlabel('f(Hz)')
 xlim([-10 10])
 ylim([-1 1])

 subplot(3,4,12)%Resultmage
imagesc(m.f,m.k,log(abs(tr1b)))
 ylabel('kx(cm-1)')
 xlabel('f(Hz)')
 xlim([-10 10])
 ylim([-1 1])