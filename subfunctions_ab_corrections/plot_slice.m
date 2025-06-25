% function plot_slices(Rxx, Rkk, Img, z_ind)
%
% plots slices of Rxx and Rkk matrix. Also plots DAS image
% Input:
% Rxx : matrix in xx basis
% Rkk : matrix in kk basis
% Img : Image parameters
% z_ind : optional. z index (depth) to plot slice. If not specified, plots
% all depth z. 
%
% Megumi Hirose (18 June 2025)
%

function plot_slice(Rxx, T0, Img, z_ind)
if nargin < 4
    z_ind = 1:Img.Nz;
end

num_kx = length(T0(1,:));
Rkk = zeros(num_kx, num_kx,Img.Nz);
% we have:  Rxx(x_in,x_out,z)
for ii = 1:Img.Nz
    Rkk(:,:,ii) = T0'*squeeze(Rxx(:,:,ii))*T0;
end

figure;
for ii=z_ind
subplot(2,2,1)
imagesc(Img.xvec*1000,Img.xvec*1000,abs(squeeze(Rxx(:,:,ii))));
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
title('Rxx')
axis image

subplot(2,2,3)
imagesc(abs(squeeze(Rkk(:,:,ii))));
title('Rkk')
colorbar
axis image
% clim([1000 6000]);
pause(0.05)
drawnow

subplot(2, 2, [2,4])
imagesc(Img.xvec*1000,Img.zvec*1000,abs(Rxx_DAS(Rxx))); % plot a 2D image of the medium sound speed
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
a = colorbar;
a.Label.String = '[m/s]'; % units for colorbar
yline(Img.zvec(ii)*1000,'w--','Linewidth',2);
axis('image'); % sets the aspect ratio so that equal tick mark increments on the x-,y- and z-axis are equal in size.
sgtitle(['Depth: ' num2str(Img.zvec(ii)*1000) ' mm,  index ' num2str(ii)]);

end