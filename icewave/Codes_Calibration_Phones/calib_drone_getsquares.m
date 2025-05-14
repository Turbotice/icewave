function [X,Y]=calib_drone_getsquares(im,BW)


Cmin = 0.2;
Cmax = 1;
Amin = 100;
Amax = 100000;

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

getframe();
