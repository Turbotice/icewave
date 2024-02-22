function movie_velocity_field(V,facq_x,facq_t,caxis_amp,fps,fig_name)

% This function create a movie of a given velocity field V
% It takes as arguments :
% - V : velocity field [nx,ny,nt]
% - facq_x : frequency for spacial sampling in box / meter
% - caxis_amp : the min and max values of the colorbar 1x2
% - fps : fps whith which the image will be displayed 
% - fig_name : string, under which the video will be saved 

video_filename = [fig_name '.avi']; % folder where the video is saved
vid = VideoWriter(video_filename);
vid.FrameRate = fps;
open(vid)

[nx,ny,nt] = size(V);
x = (1:1:nx)./facq_x;
y = (ny:-1:1)./facq_x;
t = (1:1:nt)./facq_t;
[X,Y]=meshgrid(x,y);

velocity_fig = figure;
velocity_fig.Color = [1 , 1 ,1];

for i = 1 : nt
%     V = sqrt(m.Vx(:,:,i).^2 + m.Vy(:,:,i).^2);
    
    pcolor(x,y,V(:,:,i)')
    shading interp
    title(['Time ' num2str(t(i)) ' s'],'Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    ax = gca;
    ax.FontSize = 13;
    axis image

    cbar = colorbar();
    caxis(caxis_amp)
    cbar.Label.String = '$ V \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
    cbar.FontSize = 13;
    if save_video
        T(i)=getframe(gcf);     
    end

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