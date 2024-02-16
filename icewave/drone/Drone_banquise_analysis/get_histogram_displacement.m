function get_histogram_displacement(Vx,W,fig_folder,font_size)

% This functions plots the Histogram of the pixel displacements after a PIV
% process
% This function can be used to check if the PIV settings were correct (
% window size sufficiently large, time step sufficiently big)

% It takes as arguments : 
% - Vx : the wavefield fo interest [ny,nx,nt]
% - W : the last window size used
% - fig_folder : folder where the histogram is saved

Vxmoy = nanmean(abs(Vx),3);

% first criteria : Vy > 0.1 pixels
low_bound = log10(0.1);

% second criteria : Vy < W/4
up_bound = log10(W/4);

hist_figure = figure;
hist_figure.Color = [1,1,1];
h = histogram(log10(Vxmoy(:)));
xline(low_bound,'r--')
xline(up_bound,'r--')
xlabel('$\log_{10}( \vert V_x \vert )$','Interpreter','latex','FontSize',font_size);
ylabel('$N( \vert V_x \vert )$','Interpreter','latex','FontSize',font_size);
%title('Histogram of horizontal displacements','FontSize',font_size);

% set correctly the image position for a pdf format 
set(hist_figure,'Units','Inches');
pos = get(hist_figure,'Position');
set(hist_figure,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hist_figure,'filename','-dpdf','-r0')

fig_file = [fig_folder 'histogramm_displacements'];

saveas(hist_figure, fig_file, 'fig');

end 