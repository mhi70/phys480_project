% function plot_Rxx(data_file, Rxx, Img, title_text)
%
% plots Rxx in KdB and abs
% Input: 
% data_file: file name of file containing simulated data
% Rxx : 3D matrix  Rxx(emission,detection,time)
% Returns DAS_img of Rxx matrix for plotting


function plot_Rxx(data_file, Rxx, Img, title_text)
if nargin < 3
    title_text = 'Rxx';
end
load(data_file);

t = (kgrid.t_array + pulse.length_s/2)*1e6; % time [mu s]

rx = target.dx; % target position, depth [mm]
ry = target.dy; % target position, lateral direction [mm]

img_data = Rxx_DAS(Rxx);

KdB_bmode = 20*log10(abs(img_data).^2/max(abs(img_data(:))).^2);
fig2 = figure;
abs_speed = subplot(1, 2, 1);
imagesc(Img.xvec*1000, Img.zvec*1000, KdB_bmode)
colorbar
clim([-50 0])

z_target = rx + kgrid.x_vec(1) * 1000 + abs(kgrid.x_vec(1)-array.element.depth(1))*1000 + array.element.height*1000/2;
x_target = ry + kgrid.y_vec(1) * 1000;

r = 0.45*pulse.wavelength_ref*1e3;
plot_circle_x = @(t) r*sin(t) + x_target;
plot_circle_z = @(t) r*cos(t) + z_target;
hold on;
fplot(plot_circle_x,plot_circle_z,'r:','LineWidth',1.5);
xlabel('lateral position (mm)');
ylabel('depth (mm)');
axis image;

KdB_bmode = 20*log10(abs(img_data).^2/max(abs(img_data(:))).^2);

ax_speed = subplot(1,2,2);
imagesc(Img.xvec*1000, Img.zvec*1000,abs(img_data)); % plot a 2D image of the medium sound speed
hold on;
fplot(plot_circle_x,plot_circle_z,'r:','LineWidth',1.5);
title(title_text);
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
a = colorbar;
a.Label.String = '[m/s]'; % units for colorbar
axis image
