function [xline,yline]= find_line(Xout,Yout,x1,y1,x2,y2,T)

dx = x2-x1;
dy = y2-y1;

d0 = sqrt(dx^2+dy^2);
a = sqrt(dy^2/(dx^2+dy^2));
b = -sign(dx/dy)*sqrt(dx^2/(dx^2+dy^2));
c = -a*x1-b*y1;
indices = find(abs(a*Xout+b*Yout+c)<T);

s = sqrt((Xout(indices)-x1).^2+(Yout(indices)-y1).^2);
[~,inds] = sort(s);
xline = Xout(indices(inds));
yline=Yout(indices(inds));

