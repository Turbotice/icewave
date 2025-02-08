f=rdir('calib_1/*1350.tif*')

img=imread(f{1});
filename=f{1};
filename=strcat(filename(1:end-4),'mat')
imshow(img)
%%

%filtrage en fourier

Y=fftshift(fft2(img(:,:,1),4096,4096));
size_filt=40;
Y(2049-size_filt:2049+size_filt,2049-size_filt:2049+size_filt)=zeros;

imgfilt=ifft2(ifftshift(Y));
imgfilt=imgfilt(1:size(img,1),1:size(img,2));

imagesc(imgfilt)


se = strel('disk',3)
imgfilt=imgfilt+imtophat(imgfilt,se)-imbothat(imgfilt,se);

%% Binarization
BW=imgfilt;
BW(find(imgfilt>-8))=1;
BW(find(imgfilt<=-8))=0;

%T = adaptthresh(BW);
%BW = imbinarize(BW,T);

%BW(find(imgfilt>-5))=1;
%BW(find(imgfilt<=-5))=0;
% 
SE = strel('disk',1)
% 
% %BW2=medfilt2(BW);
% BW2=imdilate(BW2,SE);
% BW2=imdilate(BW2,SE);
% 
% 
 BW2=imerode(BW,SE);
% BW2=imerode(BW2,SE);
% 
BW2 = imfill(BW2,'holes');
imshow(BW2)

%%
BW=BW2;
save(filename,'BW')

%%

% % % 
% % % BW = imfill(BW,'holes');
% % % 
% % % imshow(BW)
% % % pause
% % % L=bwlabel(BW);
% % % S=regionprops(L,'Area','Centroid','Eccentricity')
% % % 
% % % 
% % % %%
% % % %BW=BW(:,:,3);
% % % centroids = cat(1, S.Centroid);
% % % 
% % % imshow(BW)
% % % hold on
% % % plot(centroids(:,1),centroids(:,2), 'b*')
% % % hold off