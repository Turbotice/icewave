function [Tics,Average]=radialavg2(data,radial_step,x0,y0)
% This function computes the radial average of a 2D matrix
% It takes as arguments :
% - data, the 2D matrix
% - radial_step, the radial step used to average
% - x0, the shift applied on the x-axis to perform the averaging
% - y0, the shift applied on the y-axis to perform the averaging

% It returns :
% - Tics, the mean radius of each circular layer
% - Average, the value of the radial average at a given circular layer

%main axii cpecified:
l=round(size(data,1)/2);
x=(1:size(data,2))-size(data,2)/2+(l-y0);
y=(1:size(data,1))-size(data,1)/2+(l-x0);
% coordinate grid:
[X,Y]=meshgrid(x,y);
%X=X-x0;
%Y=Y-y0;

% creating circular layers
Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
Average=accumarray(Z_integer(:),data(:),[],@mean);
end