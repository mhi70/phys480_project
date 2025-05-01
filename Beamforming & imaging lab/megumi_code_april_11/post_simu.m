% process full matrix capture simulation data
% filter (Hannint + Tukey) and perform B-mode imaging
% author: Megumi Hirose
% created 9 April 2025

% load 'Kuu.mat'
% if 'Kuu.mat' not available, run 'run_simulation_full_matrix.m' first
clear variables
clear global
close all hidden
clc
load('Kuu_new.mat')

filter = true;
plot = true;

t = (kgrid.t_array + pulse.length_s/2)*1e6; % time [mu s]
Dt = max(t);
freq = (0:(kgrid.Nt-1))/Dt; % [MHz]

rx = 11; % target position, depth [mm]
ry = 12; % target position, lateral direction [mm]

%% Plot grid geometry and medium characteristics - to check final image
fig1 = figure;
ax_density = subplot(1,2,1);
imagesc(kgrid.y_vec*1000,kgrid.x_vec*1000,medium.density); % plot a 2D image of the medium density
title('density');
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
a = colorbar;
a.Label.String = 'kg/m^3'; % units for colorbar
axis image  % sets the aspect ratio so that equal tick mark increments on the x-,y- and z-axis are equal in size.

ax_speed = subplot(1,2,2);
imagesc(kgrid.y_vec*1000,kgrid.x_vec*1000,medium.sound_speed); % plot a 2D image of the medium sound speed
title('sound speed')
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
a = colorbar;
a.Label.String = '[m/s]'; % units for colorbar
axis image


%%
if filter
    for ii = 1:array.element.num % for every emitter
        K = squeeze(M(ii,:,:));
        K = K - mean(K, 2); % get rid of any constant offset
        %% Remove emission pulse (only want reflected pulse)
        for kk=1:array.element.num % loop over all detectors
            detector_x = array.element.lateral(kk); % [m]
            emitter_x = array.element.lateral(ii); % [m]
            distance = abs(detector_x - emitter_x); % depth is the same for all
            t_direct = pulse.length_s + 2*distance/medium.sound_speed_ref; % travel time straight from emitter to detector [s]
            npts_direct = round(t_direct*pulse.Fs); % pulse.Fs = samples/second ==> number of samples
            K(kk, 1:npts_direct) = 0;
        end

        %% Convert K to decibels
        K_dB = 20*log10(abs(K).^2/max(abs(K(:))).^2);

        % zpixels = (-2:0.05:6)/1000; % [m]
        % figure;
        % xpixels = (-6:0.05:6)/1000; % [m]
        % imagesc(xpixels*1000, zpixels*1000, K_dB);
        % axis image;
        % xlabel('lateral position (mm)');
        % ylabel('depth (mm)');
        % title(['Emission ' num2str(ii)]);
        % drawnow;

        %% Perform frequency filtering - Hanning window

        Img.fc = pulse.tone_burst_freq_HF;
        Img.freq_lim = [0.1, 7]; % [MHz]

        idxFreq=find(freq>=Img.freq_lim(1) & freq<=Img.freq_lim(2));
        Nfkeep = length(idxFreq); % num points to keep

        % Hanning Filter on selected frequency
        HANN=zeros(kgrid.Nt, 1);
        HANN(idxFreq)=hann(Nfkeep);

        %% Perform time domain filtering - Tukey window
        % Apply Hanning Filter to the entire dataset
        Kfft = fft(K,[],2);
        Kfft = Kfft.*HANN.';

        % Go back into the time domain
        Kfilt=ifft(Kfft,[],2);

        % Apply time-domain filter
        TimeFilt = tukeywin(kgrid.Nt,0.025);

        Kfilt_time = Kfilt.*TimeFilt.';
        M(ii, :, :) = Kfilt_time;

    end
end

%% B-mode imaging - VECTORIZED

figure;

% target position
z_target = rx - abs(kgrid.x_vec(1)-array.element.depth(1) + array.element.height/2)*1000;
x_target = ry + kgrid.y_vec(1)*1000;

% Image parameters
zpixels = (-2:0.05:6)/1000; % [m]
xpixels = (-6:0.05:6)/1000; % [m]

Nz = length(zpixels); % depth (y)
Nx = length(xpixels); % lateral (x)

img_data = zeros(Nz, Nx);

[X, Z] = meshgrid(xpixels, zpixels);  % [Nz x Nx]

for ii = 1:array.element.num
    emitter_x = array.element.lateral(ii);
    emitter_z = array.element.depth(ii); % [m]

    for xx = 1:array.element.num
        element_x = array.element.lateral(xx);
        element_z = array.element.depth(xx);

        % Get the signal for this emitter and receiver
        R = squeeze(M(ii, xx, :));  % size: [Nt x 1]

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
            imagesc(xpixels*1000, zpixels*1000, abs(img_data));
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
imagesc(xpixels*1000, zpixels*1000, KdB_bmode)
colorbar
clim([-50 0])

% target position

% z_target = rx - abs(kgrid.x_vec(1)-array.element.depth(1) + array.element.height/2)*1000;
% x_target = ry + kgrid.y_vec(1)*1000;

% % plot a circle representing the original pressure source
% source_depth_reltoarray = source_depth -source_radius;  % [mm]
% source_lateral_reltoarray = source_lateral - kgrid.y_vec(end)*1000-source_radius;  % [mm]

% rx = 11; % target position, depth [mm]
% ry = 12; % target position, lateral direction [mm]

z_target = rx + kgrid.x_vec(1) * 1000 + abs(kgrid.x_vec(1)-array.element.depth(1))*1000 + array.element.height*1000/2;
x_target = ry + kgrid.y_vec(1) * 1000;

r = 0.45*pulse.wavelength_ref*1e3;
plot_circle_x = @(t) r*sin(t) + x_target;
plot_circle_z = @(t) r*cos(t) + z_target;
hold on;
fplot(plot_circle_x,plot_circle_z,'k:','LineWidth',1.5);
xlabel('lateral position (mm)');
ylabel('depth (mm)');


