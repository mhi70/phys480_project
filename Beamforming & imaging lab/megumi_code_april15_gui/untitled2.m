% Aberration correction to US simulation
%
% Megumi Hirose (25 April 2025)

clear variables
clear global
close all hidden
clc

load('K_US_PA_aberration2.mat')
filter = true;

t = kgrid.t_array*1e6; % time [mu s]
Dt = max(t);
freq = (0:(kgrid.Nt-1))/Dt; % [MHz]

%% image params used in US B-mode images

xpixels = (-7.5:0.05:7.5)/1000; % [m] depth
ypixels = (-10.5:0.05:10.5)/1000; % [m] lateral

% target location on image
x_target = (rx + kgrid.x_vec(1) + abs(kgrid.x_vec(1) - array.element.depth(1)) + array.element.height/2)*1000; %[mm]
y_target = (ry + kgrid.y_vec(1)) * 1000; % [mm]

% for plotting target
r = 0.5*pulse.wavelength_ref*1000; % [mm]
plot_circle_x = @(t) r*sin(t) + x_target; % [mm]
plot_circle_y = @(t) r*cos(t) + y_target; % [mm]


%% Plot grid geometry and medium characteristics - to check final image
figure;
ax_density = subplot(1,2,1);
medium = medium_US;
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


%% apply filter

if filter
    for ii = 1:array.element.num % for every emitter
        K = squeeze(K_US(ii,:,:));
        K = K - mean(K, 2); % get rid of any constant offset
        % Remove emission pulse (only want reflected pulse)
        for kk=1:array.element.num % loop over all detectors
            detector_x = array.element.lateral(kk); % [m]
            emitter_x = array.element.lateral(ii); % [m]
            distance = abs(detector_x - emitter_x); % depth is the same for all
            t_direct = pulse.length_s + 2*distance/medium.sound_speed_ref; % travel time straight from emitter to detector [s]
            npts_direct = round(t_direct*pulse.Fs); % pulse.Fs = samples/second ==> number of samples
            K(kk, 1:npts_direct) = 0;
        end

        % Convert K to decibels
        K_dB = 20*log10(abs(K).^2/max(abs(K(:))).^2);

        % Perform frequency filtering - Hanning window

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
        K_US(ii, :, :) = Kfilt_time;

    end
end

% K_US = original matrix with all the data

% %% START ABERRATION CORRECTION

%% Step 1: Create Rxx
Ny = length(ypixels); % lateral (y) 421
Nx = length(xpixels); % depth (x) 301

Rxx = zeros(Ny,Ny,Nx); % input lateral,output letaral ,depth 421 x 421 x 301
DAS_img = zeros(Ny,Nx);

[X,Y] = meshgrid(xpixels,ypixels); % X, Y = 421 x 301

for ii = 1:array.element.num % for every emitter
    emitter_y = array.element.lateral(ii); % [m]
    emitter_x = array.element.depth(ii); % [m]

    DAS_ii = zeros(Ny,Nx); % DAS image at this value 421 x 301

    for xx = 1:array.element.num % for every receiver
        element_y = array.element.lateral(xx); % [m]
        element_x = array.element.depth(xx); % [m]

        % Get the signal for this emitter and receiver
        R = squeeze(K_US(ii, xx, :));  % size: [Nt x 1] 5391 x 1

        % Calculate transmit and receive distances for all pixels at once
        transmit_distance = sqrt((emitter_y - Y).^2 + (emitter_x - X).^2);  % [Ny x Nx] 421 x 301
        receive_distance = sqrt((element_y - Y).^2 + (element_x - X).^2);   % [Ny x Nx] 421 x 301

        % Compute total time delay
        time_delay = (transmit_distance + receive_distance) / medium_US.sound_speed_ref;  % [Ny x Nx] 421 x 301
        time_delay_us = time_delay * 1e6;  % Convert to microseconds

        % Interpolation: replace NaN values with 0 after interpolation
        R_interp = interp1(t, R, time_delay_us, 'linear', 0);  % [Ny x Nx] 421 x 301

        % DAS_img = DAS_img + R_interp;
        DAS_ii = DAS_ii + R_interp;
    end
    DAS_img = DAS_img + DAS_ii;
    Rxx(ii,:,:) = DAS_ii; % [Ny x Ny x Nx] 421 x 421 x 301
end


%% Recreate DAS_img with Rxx

img_data = squeeze(sum(Rxx));
KdB_bmode_US = 20*log10(abs(img_data).^2/max(abs(img_data(:))).^2);
figure;
imagesc(ypixels*1000, xpixels*1000, KdB_bmode_US)
colorbar
clim([-50 0])

hold on;
fplot(plot_circle_y,plot_circle_x,'k:','LineWidth',1.5);
xlabel('lateral position (mm)');
ylabel('depth (mm)');
title('US full matrix')
axis image

%% Plotting a slice of Rxx at depth z
z = 10;
img_data = squeeze(Rxx(:,:,z));
KdB_bmode_US = 20*log10(abs(img_data).^2/max(abs(img_data(:))).^2);
figure;
imagesc(ypixels*1000, xpixels*1000, KdB_bmode_US)
colorbar
clim([-50 0])

hold on;
fplot(plot_circle_y,plot_circle_x,'k:','LineWidth',1.5);
xlabel('lateral position (mm)');
ylabel('depth (mm)');
title('US full matrix')
axis image