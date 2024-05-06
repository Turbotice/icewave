function get_histogram_displacement(Vx,W,bool_average,font_size,filename)

% This functions plots the Histogram of the pixel displacements after a PIV
% process
% This function can be used to check if the PIV settings were correct (
% window size sufficiently large, time step sufficiently big)

% It takes as arguments : 
% - Vx : the wavefield of interest [ny,nx,nt]
% - W : the last window size used
% - bool_average : boolean to choose to average over time 
% - fig_folder : folder where the histogram is saved

% average over time, of the absolute velocity
Vxmoy = nanmean(abs(Vx),3);

% first criteria : V > 0.1 pixels
low_bound = log10(0.1);

% second criteria : V < W/4
up_bound = log10(W/4);

hist_figure = figure;
hist_figure.Color = [1,1,1];
if bool_average
    h = histogram(log10(Vxmoy(:)));
else
    h = histogram(log10(abs(Vx)));
end

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

fig_file = filename;

saveas(hist_figure, fig_file, 'fig');

end 