clear all
close all

%%
% f=rdir('Images/*.tiff')
path2tiff = 'W:/Banquise/Calibration_PMMH_2026/20251216/calib_video_4K/frame_selection/';
f = dir([path2tiff 'Binarisation/' '*.tiff']);

path2save = [ path2tiff 'Binarisation/'];
if ~isfolder(path2save) 
    mkdir(path2save)
end

%%
for i=1:length(f)
    
    disp(strcat(num2str(i),'/',num2str(length(f))))
    
    file2load = f(i).name;
    img=imread([f(i).folder '/' file2load]);
    filename=file2load;
    filename=filename(end-8:end-5);
    filename=strcat(path2save,filename,'.mat');

%filename=strcat('Binarisation/',filename(1:end-4),'mat');
%imshow(img)

    %filtrage en fourier
    padding=4096;
    Y=fftshift(fft2(img(:,:,1),padding,padding));
    size_filt=40;
    Y(padding/2+1-size_filt:padding/2+1+size_filt,padding/2+1-size_filt:padding/2+1+size_filt)=zeros;
    
    imgfilt=ifft2(ifftshift(Y));
    imgfilt=imgfilt(1:size(img,1),1:size(img,2));
    
    imagesc(imgfilt)
    
    
    se = strel('disk',3);
    imgfilt=imgfilt+imtophat(imgfilt,se)-imbothat(imgfilt,se);
    imagesc(imgfilt)

    % Binarization
    BW=imgfilt;
    BW(find(imgfilt>-8))=1;
    BW(find(imgfilt<=-8))=0;
    
    %T = adaptthresh(BW);
    %BW = imbinarize(BW,T);
    
    %BW(find(imgfilt>-5))=1;
    %BW(find(imgfilt<=-5))=0;
    % 
    SE = strel('disk',1);
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
    
    %
    BW=BW2;
    % save(filename,'BW')

end


%%
clear
f=rdir('Bin*/*.mat')

for i=1:length(f)
    load(f{i})
    
    imshow(BW)
    title(num2str(i))
    pause
   
    
    end
    
