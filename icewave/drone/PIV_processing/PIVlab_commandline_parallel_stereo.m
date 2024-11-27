% Example script how to use PIVlab from the commandline
% Just run this script to see what it does.
% You can adjust the settings in "s" and "p", specify a mask and a region of interest
function [u,v,xpix,ypix] = PIVlab_commandline_parallel_stereo(directory,reference,N,i0,Dt,b,s,p)

% Create list of images inside specified directory
suffix='*.tiff'; %*.bmp or *.tif or *.jpg or *.tiff or *.jpeg
filesep = '/';
direc = dir ( [directory filesep suffix] ); filenames={};

[filenames{1:length(direc),1}] = deal(direc.name);
filenames = sortrows(filenames); %sort all image files
amount = length(filenames);

disp(amount)


% PIV analysis loop


if mod(amount,2) == 1 %Uneven number of images?
    disp('Image folder should contain an even number of images.')
    %remove last image from list
    amount=amount-1;
    filenames(size(filenames,1))=[];
end
disp(['Found ' num2str(amount) ' images (' num2str(amount-1) ' image pairs).'])
x=cell(amount-1,1);
y=x;
u=x;
v=x;
typevector=x; %typevector will be 1 for regular vectors, 0 for masked areas
typevector0=1;
% counter=0;
% PIV analysis loop:
%N = 3000;%amount-1;%floor(amount/2);

filenames
if ~strcmp(reference,'')
    filename1 = reference;
    %filename1 = filenames{1};
    filename2 = filenames{1};
else
    filename1 = filenames{1};
    filename2 = filenames{1+Dt};
end

%image1=imread(fullfile(directory, filename1)); % read images
%image2=imread(fullfile(directory, filename2));

nr_of_cores = 20;

[x0, y0, u0, v0, typevector,correlation_map] = ...
			piv_analysis(directory, filename1, filename2,p,s,nr_of_cores,false);

if N==0
    N = amount;
end

Nmax = floor((N-Dt-1-i0)/b+1);

U = zeros([Nmax size(u0)]);
V = zeros([Nmax size(v0)]);

disp(['Number of images to process : ' num2str(Nmax)])

parfor (i=1:Nmax,nr_of_cores)
% parfor = for loop processed on parallel cores 
    if ~strcmp(reference,'')
        filename1 = reference;
        %filename1 = filenames{1};
        filename2 = filenames{i0+(i-1)*b+1};
    else
        filename1 = filenames{i0+(i-1)*b+1};
        filename2 = filenames{i0+(i-1)*b+1+Dt};
    end
    image1=imread(fullfile(directory, filename1)); % read images
    image2=imread(fullfile(directory, filename2));
    
    image1 = PIVlab_preproc (image1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2},p{9,2},p{10,2}); %preprocess images
    image2 = PIVlab_preproc (image2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2},p{9,2},p{10,2});
    [~, ~, u0, v0, ~] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2},s{11,2},s{12,2},0,0,0,0);
    
    U(i,:,:) = u0;
    V(i,:,:) = v0;
    
    disp([int2str((i-i0)/Nmax*100) ' %']);
    
    % Graphical output (disable to improve speed)
    %%{
    %     imagesc(double(image1)+double(image2));colormap('gray');
    %     hold on
    %     quiver(x{counter},y{counter},u{counter},v{counter},'g','AutoScaleFactor', 1.5);
    %     hold off;
    %     axis image;
    %     title(filenames{i},'interpreter','none')
    %     set(gca,'xtick',[],'ytick',[])
    %     drawnow;
    %%}
end

% counter = 0; % ?
% for i=1:1:Nmax
%     counter = counter+1;
%     x{counter} = x0;
%     y{counter} = y0;
%     u{counter} = squeeze(U(counter,:,:));
%     v{counter} = squeeze(V(counter,:,:));
%     typevector{counter} = typevector0;
% end

for i=1:1:Nmax
%     counter = counter+1;
    x{i} = x0;
    y{i} = y0;
    u{i} = squeeze(U(i,:,:));
    v{i} = squeeze(V(i,:,:));
    %typevector{i} = typevector0;
end

% create array of box position 
xpix = x0(1,:); % x-pixel position of center of each box
ypix = y0(:,1); % y-pixel position of center of each box 

%% PIV postprocessing loop
% Settings
if false
    W = 16;
    umin = -W/2; % minimum allowed u velocity, adjust to your data
    umax = W/2; % maximum allowed u velocity, adjust to your data
    vmin = -W/2; % minimum allowed v velocity, adjust to your data
    vmax = W/2; % maximum allowed v velocity, adjust to your data
    stdthresh=7; % threshold for standard deviation check
    epsilon=0.15; % epsilon for normalized median test
    thresh=3; % threshold for normalized median test
    u_filt=cell(N,1);
    v_filt=u_filt;
    typevector_filt=u_filt;
    for PIVresult=1:size(x,1)
        u_filtered=u{PIVresult,1};
        v_filtered=v{PIVresult,1};
        typevector_filtered=typevector{PIVresult,1};
        %vellimit check
        u_filtered(u_filtered<umin)=NaN;
        u_filtered(u_filtered>umax)=NaN;
        v_filtered(v_filtered<vmin)=NaN;
        v_filtered(v_filtered>vmax)=NaN;
        % stddev check
        meanu=nanmean(nanmean(u_filtered));
        meanv=nanmean(nanmean(v_filtered));
        std2u=nanstd(reshape(u_filtered,size(u_filtered,1)*size(u_filtered,2),1));
        std2v=nanstd(reshape(v_filtered,size(v_filtered,1)*size(v_filtered,2),1));
        minvalu=meanu-stdthresh*std2u;
        maxvalu=meanu+stdthresh*std2u;
        minvalv=meanv-stdthresh*std2v;
        maxvalv=meanv+stdthresh*std2v;
        u_filtered(u_filtered<minvalu)=NaN;
        u_filtered(u_filtered>maxvalu)=NaN;
        v_filtered(v_filtered<minvalv)=NaN;
        v_filtered(v_filtered>maxvalv)=NaN;
        % normalized median check
        %Westerweel & Scarano (2005): Universal Outlier detection for PIV data
        [J,I]=size(u_filtered);
        medianres=zeros(J,I);
        normfluct=zeros(J,I,2);
        b=1;
        for c=1:2
            if c==1; velcomp=u_filtered;else;velcomp=v_filtered;end %#ok<*NOSEM>
            for i=1+b:I-b
                for j=1+b:J-b
                    neigh=velcomp(j-b:j+b,i-b:i+b);
                    neighcol=neigh(:);
                    neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
                    med=median(neighcol2);
                    fluct=velcomp(j,i)-med;
                    res=neighcol2-med;
                    medianres=median(abs(res));
                    normfluct(j,i,c)=abs(fluct/(medianres+epsilon));
                end
            end
        end
        info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
        u_filtered(info1==1)=NaN;
        v_filtered(info1==1)=NaN;
        typevector_filtered(isnan(u_filtered))=2;
        typevector_filtered(isnan(v_filtered))=2;
        typevector_filtered(typevector{PIVresult,1}==0)=0; %restores typevector for mask
        
        %Interpolate missing data
        u_filtered=inpaint_nans(u_filtered,4);
        v_filtered=inpaint_nans(v_filtered,4);
        
        u_filt{PIVresult,1}=u_filtered;
        v_filt{PIVresult,1}=v_filtered;
        typevector_filt{PIVresult,1}=typevector_filtered;
        disp('Entered in the final loop');
    end
end

%clearvars -except p s x y u v typevector directory filenames u_filt v_filt typevector_filt
disp('DONE.')
end
