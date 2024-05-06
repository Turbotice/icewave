function [profile_fig] = plot_located_profile(V,i_x,i_y,fx,fps,left_bool)

% This function computes the time profile of the velocity at a given
% position i_x,i_y

% - V : velocity field [nx,ny,nt]
% - i_x : index along the x-axis (dimension 1)
% - i_y : index along the y-axis (dimension 2)
% - fx : spacial scale 
% - fps : time scale 
% - left_bool : boolean to choose how to orient x-axis
% if left_bool = 1 -> x-axis positive from left to right
% if left_bool = 0 -> x-axis positive from right to left

[nx,ny,nt] = size(V);
if left_bool 
    x = (1:1:nx);
else 
    x = (nx:-1:1);
end 
y = (1:1:ny);
t = (1:1:nt);
[X,Y] = meshgrid(x,y);

profile_fig = figure;
profile_fig.Color = [1, 1, 1];

subplot(2,1,1)
plot(t/fps,squeeze(V(i_x,i_y,:)));
xlim([0 nt/fps]);
xlabel('$t \: \rm (s)$','Interpreter','latex');
ylabel('$V \: \rm (m/s)$','Interpreter','latex');

ax = gca;
ax.FontSize = 13;

subplot(2,1,2)
pcolor(x/fx,y/fx,V(:,:,1)')
shading interp
hold on 
plot(x(i_x)/fx,y(i_y)/fx,'ro')
xlabel('$x \: \rm (m)$','Interpreter','latex');
ylabel('$y \: \rm (m)$','Interpreter','latex');

if ~left_bool
    set(gca,'XDir','reverse')
end 

ax = gca;
ax.FontSize = 13;

end 