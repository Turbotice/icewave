
M={};

folder = 'Phone_01/';

imlist = dir([folder '*.jpg']);
disp(imlist)

%%
c=0;
keep=[];
for i=1:length(imlist)
    imname = imlist(i).name;
    im = imread(fullfile(folder,imname));

    imshow(im)
    getframe()
    s = input('Keep image ?')
    if length(s)==0
        disp('good !')
        c=c+1;
        keep(c)=i;
    end
end

%%
keep_04 = keep;
%%




folder = 'Phone_01/';

imlist = dir([folder '*.jpg']);
disp(imlist)

for i=1:length(imlist)
    if sum(keep_01==i)>0
        imname = imlist(i).name;
        filename = fullfile(folder,imname);
        im = imread(filename);
        filename=strcat(filename(1:end-4));%,'mat');
        BW = binarize_smart(im);

        save(filename,'BW')
    end
end




%%

%folder = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/Calibration_intrinseque/mesange/selection/calib_1/';%
%folder = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/Calibration_intrinseque/bernache/calib_video/Binarisation/';

folder = 'Phone_01/';

bwlist = dir([folder '*.mat']);
%imlist = dir([folder '*.jpg']);

figure(1)
M={};
for i=1:length(bwlist)
    name = bwlist(i).name;
    cel = split(name,'.');
    imname = [cel{1} '.jpg'];
    %imname = imlist(i).name;

    disp(name)
    disp(imname)
    load(fullfile(bwlist(i).folder,name))
    im = imread(fullfile(bwlist(i).folder,imname));

    imshow(BW)
    %title(name)

    m.name = name;
    %m.corners = calib_drone_getcorners(im);

    getframe();

    s = input("");

    
    M(i).m=m;
    M(i).m.imname = imname;
end
%
%Mref = M;

%%
for i=1:length(M)
    name = M(i).m.name;
    imname = M(i).m.imname;
    load(fullfile(folder,name))
    im = imread(fullfile(imlist(i).folder,imname));

    figure(5)
    imshow(BW)
    hold on

    [X,Y] = calib_drone_getsquares(im,BW);
    %    [X,Y]=calib_drone_getgrid(im,BW);

    M(i).m.X=X;
    M(i).m.Y=Y;
    title(name)
    %pause(2)
end

%%

for i=1:19%length(M)-2
    name = M(i).m.name;
    imname = M(i).m.imname;
    load(fullfile(folder,name))
    im = imread(fullfile(imlist(i).folder,imname));

    figure(5)
    imshow(BW)
    hold on


    X = M(i).m.X;
    Y = M(i).m.Y;

    Msub = calib_drone_getsubgrid(im,X,Y,replace(imname,'.','_'));
    %    [X,Y]=calib_drone_getgrid(im,BW);
    M(i).m.Msub=Msub;

    title(name)
    pause(2)
end

%%

%%

matlist = dir('transient_IMG_*.mat');

Mp = {};
c=0;
for i=1:19%length(Mref)-2
    load(matlist(i).name)

    for k=1:length(M)
        n = length(M(k).m.Xsub);
        disp(n)
        if n==35
            c=c+1;
            Mp(c).m=Mref(i).m;
            Mp(c).m.X = M(k).m.Xsub;
            Mp(c).m.Y = M(k).m.Ysub;
        end
    end
    disp(c)
end


%%


%%
filename = 'matrix_points_phone_00.mat';
save(filename,'Mp')

%%
nx = 7;
ny = 5;

figure(39)
ind = {};
for i=1:length(Mp)
    m = Mp(i).m;
    Xout = m.X;
    Yout = m.Y;

    xc = median(Xout);
    yc = median(Yout);

    %     disp(i)
    hold off
    plot(Xout,Yout,'x')
    hold on
    plot(xc,yc,'o')
    d = sqrt((Xout-xc).^2+(Yout-yc).^2);

    [~,ind{1}] = max(Xout+Yout-xc-yc);
    [~,ind{2}] = max(Xout-Yout-xc+yc);
    [~,ind{3}] = max(-Xout+Yout+xc-yc);
    [~,ind{4}] = max(-Xout-Yout+xc+yc);

    Mp(i).m.center=[xc,yc];

    for k=1:4
        corners{k}=[Xout(ind{k}) Yout(ind{k})];
    end

    c=0;
    distance=[];

    for k=1:4
        for j=1:4
            if and(k~=j,k>j)
                c=c+1;
                a = sum(abs(corners{k}-corners{j}));
                distance(c)=a;
                ind1(c)=k;
                ind2(c)=j;
            end
        end
    end


    if min(distance)==0
        for k=1:4
            plot(corners{k}(1),corners{k}(2),'ks')
                    text(corners{k}(1) + 2,corners{k}(2) + 2,num2str(k))

        end
        title(i)
%         if sum(distance==0)>0
%             [xa,ya] = getpts(gcf);
%             indices = find(distance==0);
%             j = ind2(indices(1));
% 
%             d = sqrt((Xout-xa).^2+(Yout-ya).^2);
% 
%             [val,inda]=min(d);
%             m.corners{j}=[Xout(inda),Yout(inda)];
%             plot(m.corners{j}(1),m.corners{j}(2),'rs')
%         end
        getframe();
    end
    %Mp(i).m.corners = m.corners;

%

    [~,indices]=sort(-distance);
    disp(ind1(indices))
    disp(ind2(indices))

    %     i1=ind{ind1(indices(1))};
    %     i2=ind{ind1(indices(2))};
    %     i3=ind{ind2(indices(2))};
    %     i4=ind{ind2(indices(1))};

    i1=ind1(indices(1));
    i2=ind1(indices(2));
    i3=ind2(indices(2));
    i4=ind2(indices(1));


    %     disp(i1)
    %     disp(i2)
    %     disp(i3)
    %     disp(i4)
    d1 = abs(sum(corners{i1}-corners{i2}));
    d2 = abs(sum(corners{i1}-corners{i3}));

    if d1<d2
        ib=i2;
        i2=i3;
        i3=ib;
    end

    corners{1}(1)=Xout(ind{i1});
    corners{2}(1)=Xout(ind{i2});
    corners{3}(1)=Xout(ind{i3});
    corners{4}(1)=Xout(ind{i4});
    corners{1}(2)=Yout(ind{i1});
    corners{2}(2)=Yout(ind{i2});
    corners{3}(2)=Yout(ind{i3});
    corners{4}(2)=Yout(ind{i4});

    Mp(i).m.corners = corners;
    for k=1:4
        %              plot(Xout(ind{k}),Yout(ind{k}),'ks')
        plot(corners{k}(1),corners{k}(2),'ks')
        text(corners{k}(1) + 2,corners{k}(2) + 2,num2str(k))
    end
    getframe()
    %pause()
end

%%
nx = 7;
ny = 5;

lines = [[1 3];[2 4]];
lines2 = [[1 2];[3 4]];

bounds = [1 ny];

%%

corners{k}

%%


Xtab = [];
Ytab = [];
for i=1:length(Mp)


    disp(i)
    m = Mp(i).m;
    Xout = m.X;
    Yout = m.Y;

    figure(10)
    hold off
    plot(Xout,Yout,'x')

    title([num2str(i) ', ' m.name])
    hold on

    for k=1:2
        line = lines(k,:);
        x1 = m.corners{line(1)}(1);
        y1 = m.corners{line(1)}(2);

        x2 = m.corners{line(2)}(1);
        y2 = m.corners{line(2)}(2);

        T = 20;
        [xline,yline] = find_line(Xout,Yout,x1,y1,x2,y2,T);
        size(xline)
        if length(xline)==7
            Xtab(bounds(k),:) = xline;
            Ytab(bounds(k),:) = yline;
            plot(xline,yline,'ro')
        else
            for k=1:2
                line = lines2(k,:);

                x1 = m.corners{line(1)}(1);
                y1 = m.corners{line(1)}(2);

                x2 = m.corners{line(2)}(1);
                y2 = m.corners{line(2)}(2);

                T = 20;
                [xline,yline] = find_line(Xout,Yout,x1,y1,x2,y2,T);
                %size(Xtab)
                [nlinex,nliney] = size(xline);
                if nliney==7
                    Xtab(bounds(k),:) = xline;
                    Ytab(bounds(k),:) = yline;
                    plot(xline,yline,'ro')
                else
                    disp(size(xline))
                    disp('sorry')
                end
            end

        end
        getframe();
        pause(0.5)
    end

    
    for p=1:nx
        x1 = Xtab(bounds(1),p);
        y1 = Ytab(bounds(1),p);

        x2 = Xtab(bounds(2),p);
        y2 = Ytab(bounds(2),p);

        T = 25;
        [xline,yline] = find_line(Xout,Yout,x1,y1,x2,y2,T);
        disp(size(xline))
        if length(xline)==5
            Xtab(:,p) = xline;
            Ytab(:,p) = yline;
        end
        plot(xline,yline,'ro')

    end

    M(i).m.Xtab = Xtab;
    M(i).m.Ytab = Ytab;
end
%% check validity


figure(1)
c=0;
valid = [];
for i=1:length(M)
    hold off
    for p=1:nx
        plot(M(i).m.Xtab(:,p),M(i).m.Ytab(:,p),'ro')
        hold on
    end
    plot(M(i).m.Xtab(:,1),M(i).m.Ytab(:,1),'r*')
    plot(M(i).m.Xtab(:,7),M(i).m.Ytab(:,7),'r*')

    title(i)
    s=input('Valid ?')
    if length(s)==0
        c=c+1;
        valid(c)=i;
    end
end

%%


c=0;

imagePoints=[];
for i=1:length(M)
    if sum(valid==i)>0%mod(i,29)>25
        c=c+1;
        imagePoints(:,1,c) = reshape(M(i).m.Xtab,[nx*ny 1]);
        imagePoints(:,2,c) = reshape(M(i).m.Ytab,[nx*ny 1]);
    end
end
c

%%
boardSize = [6,8];
squareSizeInMM = 3955 / 40;
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
%imageSize = [3840,2160];
imageSize = [4080,3072];



%%
%étape de calibration à partir d'une référence et des images prises
stereoParams = estimateCameraParameters(imagePoints,worldPoints,ImageSize=imageSize);

%%
showReprojectionErrors(stereoParams)

figure
showExtrinsics(stereoParams)

save('stereoParams_Phone00.mat', 'stereoParams');

