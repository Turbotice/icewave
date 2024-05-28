% %% Computation of vertical displacement # test 1 
% 
% for i0 = 1:1:size(Vy,3)
%     
%     vx = Vx(:,:,i0);
%     vy = Vy(:,:,i0);
%     
%     a = H.*vx./(f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0))...  
%     - H*cos(alpha_0)*(x_pix - x_0).*vy./((f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0)).^2);
% 
%     b = H.*vy/sin(alpha_0)./(f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0)) ...
%         - H*(y_pix - y_0).*vy/tan(alpha_0)./((f*sin(alpha_0) + (y_pix - y_0).*cos(alpha_0)).^2);
%     
% %     Vz = H*( 1 - (tanalpha_0*f - (y_prime - y_0))./(f + tanalpha_0*(y_prime - y_0)) .* (Y./H + 1/tanalpha_0));
%     
%     Vz = H.*sign(b).*sqrt(a.^2 + b.^2)./(sqrt(X.^2 + (H/tan(alpha_0) - Y).^2)); 
%     dz(:,:,i0) = Vz.*fps/Dt;
%     
%     disp(i0)
%     
%     surf(X,Y,squeeze(dz(:,:,i0)))
%     shading interp
%     colorbar()
%     view(2)
%     title(num2str(i0))
%     pause(0.01)
% 
%     
% end
% 
% %% Computation of vertical displacement # test 2 
% 
% fig_folder = base;
% vid_name = 'Vertical_displacement_bernache';
% 
% save_video = 0;
% 
% if save_video
%     video_filename = [fig_folder vid_name '.avi']; % folder where the video is saved
%     vid = VideoWriter(video_filename);
%     vid.FrameRate = fps;
%     open(vid)
% end
% 
% for i0 = 1:1:size(Vy,3)
%     
%     vx = Vx(:,:,i0);
%     vy = Vy(:,:,i0);
%     
%     y_prime = y_pix + vy;
%     x_prime = x_pix + vx;
%     
%     X_prime = (x_prime - x_0)*H./(f*sin(alpha_0) + (y_prime - y_0).*cos(alpha_0));
%     Y_prime = (y_prime - y_0)*H/sin(alpha_0)./(f*sin(alpha_0) + (y_prime - y_0).*cos(alpha_0));
% %     Vz = H*( 1 - (tanalpha_0*f - (y_prime - y_0))./(f + tanalpha_0*(y_prime - y_0)) .* (Y./H + 1/tanalpha_0));
%     
%     X_prime = -X_prime;
%     Y_prime = -Y_prime;
% 
%     Vz = H.*sign(Y_prime - Y).*sqrt((X_prime - X).^2 + (Y_prime - Y).^2)./(sqrt(X_prime.^2 + (H/tan(alpha_0) + Y_prime).^2)); 
%     dz_2(:,:,i0) = Vz.*fps/Dt;
%     
%     disp(i0)
%     
%     surf(X,Y,squeeze(dz_2(:,:,i0)))
%     shading interp
%     cbar = colorbar();
%     cbar.Label.String = '$V_z \: \rm (m/s)$';
%     cbar.Label.Interpreter = 'latex';
%     caxis([-2.5 2.5])
%     view(2)
%     title(num2str(i0))
% 
%     if save_video
%         T(i0)=getframe(gcf);
%     else 
%         pause(0.01)
%     end 
%     
% end
% 
% if save_video
%     all_valid = true;
%     flen = length(T);
%     for K = 1 : flen
%       if isempty(T(K).cdata)
%         all_valid = false;
%         fprintf('Empty frame occurred at frame #%d of %d\n', K, flen);
%       end
%     end
%     if ~all_valid
%        error('Did not write movie because of empty frames')
%     end
% 
%     writeVideo(vid,T)
%     close(vid)   
% end 