% process full matrix capture simulation data
% filter (Hanning + Tukey) and perform B-mode imaging
% author: Megumi Hirose
% created 9 April 2025

function plot_Ruu(data_file, Img)
load(data_file);
plot = 0;
t = (kgrid.t_array + pulse.length_s/2)*1e6; % time [mu s]

rx = target.dx; % target position, depth [mm]
ry = target.dy; % target position, lateral direction [mm]

%% Plot grid geometry and medium characteristics - to check final image
fig1 = figure;
ax_density = subplot(1,2,1);
imagesc(kgrid.y_vec*1000,kgrid.x_vec*1000,medium_US.density); % plot a 2D image of the medium density
title('density');
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
a = colorbar;
a.Label.String = 'kg/m^3'; % units for colorbar
axis image  % sets the aspect ratio so that equal tick mark increments on the x-,y- and z-axis are equal in size.

ax_speed = subplot(1,2,2);
imagesc(kgrid.y_vec*1000,kgrid.x_vec*1000,medium_US.sound_speed); % plot a 2D image of the medium sound speed
title('sound speed')
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
a = colorbar;
a.Label.String = '[m/s]'; % units for colorbar
axis image


%% apply time and frequency filter
K_US = tukey_hann_filter(data_file);

%% B-mode imaging - VECTORIZED

% Image parameters

img_data = zeros(Img.Nz, Img.Nx);

[X, Z] = meshgrid(Img.xvec, Img.zvec);  % [Nz x Nx]

for ii = 1:array.element.num
    emitter_x = array.element.lateral(ii);
    emitter_z = array.element.depth(ii); % [m]

    for xx = 1:array.element.num
        element_x = array.element.lateral(xx);
        element_z = array.element.depth(xx);

        % Get the signal for this emitter and receiver
        R = squeeze(K_US(ii, xx, :));  % size: [Nt x 1]

        % Calculate transmit and receive distances for all pixels at once
        transmit_distance = sqrt((emitter_x - X).^2 + (emitter_z - Z).^2);  % [Nz x Nx]
        receive_distance = sqrt((element_x - X).^2 + (element_z - Z).^2);   % [Nz x Nx]

        % Compute total time delay
        time_delay = (transmit_distance + receive_distance) / medium.sound_speed_ref;  % [Nz x Nx]
        time_delay_us = time_delay * 1e6;  % Convert to microseconds

        % Interpolation: replace NaN values with 0 after interpolation
        R_interp = interp1(t, R, time_delay_us, 'linear', 0);  % [Nz x Nx]

        % Accumulate results
        img_data = img_data + R_interp;

        % --- Plotting ---
        if plot
            imagesc(Img.xvec*1000, Img.zvec*1000, abs(img_data));
            axis image;
            xlabel('lateral position (mm)');
            ylabel('depth (mm)');
            title(['Emission ' num2str(ii)]);
            drawnow;
        end
    end
end


%%
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

ax_speed = subplot(1,2,2);
imagesc(Img.xvec*1000, Img.zvec*1000,abs(img_data)); % plot a 2D image of the medium sound speed
hold on;
fplot(plot_circle_x,plot_circle_z,'r:','LineWidth',1.5);
title('Ruu')
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
a = colorbar;
a.Label.String = '[m/s]'; % units for colorbar
axis image

end
