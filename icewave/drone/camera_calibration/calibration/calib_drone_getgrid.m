function [X,Y]=calib_drone_getgrid(im,BW,corners)


Cmin = 0.4;
Cmax = 1;
Amin = 200;
Amax = 8000;

%Detect the grid points
L = bwlabel(BW,4);
S = regionprops(L,'Area','Circularity','Centroid','Perimeter');

%use double criteria on Area and circularity : better to detect less points
%than to detect bad points
%beware of the camera distance that will change area of squares

loglog([S.Circularity],[S.Area],'o')
Circ = [S.Circularity];
Area = [S.Area];

indices= find(and(and(Area>Amin,Area<Amax),and(Circ>Cmin,Circ<Cmax)));

figure(6)
hold off
imshow(im)
hold on
for i = 1:length(indices)
    pos = S(indices(i)).Centroid;
    plot(pos(1),pos(2),'ro','markersize',4,'MarkerFaceColor','r')
    X(i)=pos(1);
    Y(i)=pos(2);
    hold on
end


for j=1:length(corners)
    plot(corners{j}(1),corners{j}(2),'ks','markersize',10,'MarkerFaceColor','k')
end
getframe();

% click on image to select points 

% select only points in the rectangl defined by tape
xv = [];
yv = [];
for j = 1 : length(corners)
    xv = cat(2,xv,corners{j}(1));
    yv = cat(2,yv,corners{j}(2));
end 
xv([3 4]) = xv([4 3]);
yv([3 4]) = yv([4 3]);
in = inpolygon(X,Y,xv,yv);

figure(6)
hold off
imshow(im)
hold on
plot(X(in),Y(in),'ro','markersize',3,'MarkerFaceColor','r')
getframe();

%[vx,vy] = voronoi(X(in),Y(in));
%voronoi(X(in),Y(in))

X = X(in);
Y = Y(in);

