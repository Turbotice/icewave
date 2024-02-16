function plot_located_profile(V,i_x,i_y,fx,fps,fig_folder)

% This function computes the time profile of the velocity at a given
% position i_x,i_y

% - V : velocity field [nx,ny,nt]
% - i_x : index along the x-axis (dimension 1)
% - i_y : index along the y-axis (dimension 2)
% - fx : spacial scale 
% - fps : time scale 

[nx,ny,nt] = size(V);
x = (1:1:nx);
y = (1:1:ny);
t = (1:1:nt);
[X,Y] = meshgrid(x,y);

profile_fig = figure;
profile_fig.Color = [1, 1, 1];

subplot(2,1,1)
plot(t/fps,squeeze(V(i_x,i_y,:)));
xlim([0 nt/fps]);
xlabel('$t \: \rm (s)$','Interpreter','latex');
ylabel('$V \: \rm (pix)$','Interpreter','latex');

ax = gca;
ax.FontSize = 13;

subplot(2,1,2)
imagesc(x/fx,y/fx,V(:,:,1)')
hold on 
plot(x(i_x)/fx,y(i_y)/fx,'ro')
xlabel('$x \: \rm (m)$','Interpreter','latex');
ylabel('$y \: \rm (m)$','Interpreter','latex');

ax = gca;
ax.FontSize = 13;

fig_file_name = [fig_folder 'Time_profile_ix_' num2str(i_x) '_iy_' num2str(i_y)];
saveas(profile_fig,fig_file_name,'fig')

end 