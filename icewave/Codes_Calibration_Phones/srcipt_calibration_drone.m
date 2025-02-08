
M={};

%%

%folder = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/Calibration_intrinseque/mesange/selection/calib_1/';%

folder = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/Calibration_intrinseque/bernache/calib_video/Binarisation/';


bwlist = dir([folder '*.mat']);
imlist = dir([folder '*.jpg']);

figure(1)
for i=1:length(bwlist)-2
    name = bwlist(i).name;
    imname = imlist(i).name;

    disp(name)
    disp(imname)
    load(fullfile(bwlist(i).folder,name))
    im = imread(fullfile(imlist(i).folder,imname));

    %imshow(BW)
    %title(name)

    m.name = name;
    %m.corners = calib_drone_getcorners(im);

    %    getframe()
    %pause
    M(i).m=m;
end
%

for i=1:length(M)-2
    imname = imlist(i).name;
    M(i).m.imname = imname;
end

Mref = M;

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
    pause(2)
end

%%

for i=9:length(M)-2
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

    title(name)
    pause(2)
end

%%
Mp = {};

%%

matlist = dir('transient_im_*.mat');

c=0;
for i=1:length(Mref)-2
    load(matlist(i).name)

    for k=1:length(M)
        n = length(M(k).m.Xsub);
        if n==63
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
filename = 'matrix_points_bernache.mat';
save(filename,'Mp')

%%
nx = 9;
ny = 7;

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


    if i==107
        for k=1:4
            corners{k}=[Xout(ind{k}) Yout(ind{k})];
        end
    else
        corners{k} = Mp(i).m.corners{k};
    end

    c=0;
    distance=[];

    i
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
nx = 9;
ny = 7;

lines = [[1 3];[2 4]];
lines2 = [[1 2];[3 4]];

bounds = [1 ny];

%%

corners{k}

%%

for i=1:length(Mp)
    i
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

        T = 10;
        [xline,yline] = find_line(Xout,Yout,x1,y1,x2,y2,T);

        if length(xline)==9
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

                T = 10;
                [xline,yline] = find_line(Xout,Yout,x1,y1,x2,y2,T);
            end

            Xtab(bounds(k),:) = xline;
            Ytab(bounds(k),:) = yline;
        end
        getframe();
        pause(0.5)
    end
end

%%
%
for p=1:nx

    x1 = Xtab(bounds(1),p);
    y1 = Ytab(bounds(1),p);

    x2 = Xtab(bounds(2),p);
    y2 = Ytab(bounds(2),p);

    T = 10;
    [xline,yline] = find_line(Xout,Yout,x1,y1,x2,y2,T);

    Xtab(:,p) = xline;
    Ytab(:,p) = yline;
end

M(i).m.Xtab = Xtab;
M(i).m.Ytab = Ytab;

%%

c=0;

imagePoints=[];
for i=1:length(M)
    if mod(i,29)>25
        c=c+1;
        imagePoints(:,1,c) = reshape(M(i).m.Xtab,[nx*ny 1]);
        imagePoints(:,2,c) = reshape(M(i).m.Ytab,[nx*ny 1]);
    end
end

%%
boardSize = [8,10];
squareSizeInMM = 3955 / 40;
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
imageSize = [3840,2160];


%%
%étape de calibration à partir d'une référence et des images prises
stereoParams = estimateCameraParameters(imagePoints,worldPoints,ImageSize=imageSize);

%%
showReprojectionErrors(stereoParams)

figure
showExtrinsics(stereoParams)

save('stereoParams.mat', 'stereoParams')