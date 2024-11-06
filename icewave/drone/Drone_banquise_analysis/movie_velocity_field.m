function movie_velocity_field(V,facq_x,facq_t,step,caxis_amp,fps,fig_name)

% This function create a movie of a given velocity field V
% It takes as arguments :
% - V : velocity field [nx,ny,nt]
% - facq_x : frequency for spacial sampling in box / meter
% - facq_t : acquisition frequency for time 
% - step : step with which we go through time array
% - caxis_amp : the min and max values of the colorbar 1x2
% - fps : fps whith which the image will be displayed 
% - fig_name : string, under which the video will be saved 

video_filename = [fig_name '.avi']; % folder where the video is saved
vid = VideoWriter(video_filename);
vid.FrameRate = fps;
vid.Quality = 90; 
open(vid)

[nx,ny,nt] = size(V);

x = (1:1:nx)./facq_x;
y = (ny:-1:1)./facq_x;
t = (1:1:nt)./facq_t;

[X,Y]=meshgrid(x,y);

velocity_fig = figure;
velocity_fig.Color = [1 , 1 ,1];

relevant_idx = (1:step:nt);
for i0 = 1 : length(relevant_idx)
    i = relevant_idx(i0);
%     V = sqrt(m.Vx(:,:,i).^2 + m.Vy(:,:,i).^2);
    
    pcolor(x,y,V(:,:,i)')
    shading interp
    title(['Time ' num2str(t(i)) ' s'],'Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    ax = gca;
    ax.FontSize = 13;
    axis image

%     colormap(redblue)
    cbar = colorbar();
    caxis(caxis_amp)
    cbar.Label.String = '$ V_x \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
    
    % set figure resolution 
    set(gcf, 'Position', [100, 100, 1280, 720]);
    ax = gca;
    ax.FontSize = 16;
    cbar.FontSize = 16;
%     set_Papermode(gcf)

    T(i0)=getframe(gcf);     


end 


all_valid = true;
flen = length(T);
for K = 1 : flen
  if isempty(T(K).cdata)
    all_valid = false;
    fprintf('Empty frame occurred at frame #%d of %d\n', K, flen);
  end
end
if ~all_valid
   error('Did not write movie because of empty frames')
end

writeVideo(vid,T)
close(vid)    




end 