% plot PA data from running 'simu_kwave_convention.m'
% Megumi Hirose (12 april 2025)
% load 'K_US_PA.mat'
% if 'K_US_PA.mat' not available, run 'US_PA_simu' first

clear variables
clear global
close all hidden
clc

load('K_US_PA_aberration2.mat')

t = kgrid.t_array*1e6; % time [mu s]
Dt = max(t);
freq = (0:(kgrid.Nt-1))/Dt; % [MHz]

%% image params used in both PA and US B-mode images

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


%% plot US image

% for every emitter
K = squeeze(K_US(20,1,:));
K = K - mean(K,2); % get rid of any constant offset

Fs = length(t);
% Continuous Wavelet Transform
figure;
cwt(K,Fs) ; % This uses default Morse wavelet
ax = gca ;
ax.FontSize =25;
title('CWT of Transient Signal ') 
%%
filter = true;
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

% B-mode imaging - VECTORIZED
plot = false;

% Image parameters

Nx = length(xpixels); % depth (x)
Ny = length(ypixels); % lateral (y)

img_data = zeros(Nx, Ny);

[Y,X] = meshgrid(ypixels, xpixels);  % [Nx x Ny]

for ii = 1:array.element.num
    emitter_y = array.element.lateral(ii);
    emitter_x = array.element.depth(ii); % [m]

    for xx = 1:array.element.num
        element_y = array.element.lateral(xx);
        element_x = array.element.depth(xx);

        % Get the signal for this emitter and receiver
        R = squeeze(K_US(ii, xx, :));  % size: [Nt x 1]

        % Calculate transmit and receive distances for all pixels at once
        transmit_distance = sqrt((emitter_y - Y).^2 + (emitter_x - X).^2);  % [Nz x Nx]
        receive_distance = sqrt((element_y - Y).^2 + (element_x - X).^2);   % [Nz x Nx]

        % Compute total time delay
        time_delay = (transmit_distance + receive_distance) / medium_US.sound_speed_ref;  % [Nz x Nx]
        time_delay_us = time_delay * 1e6;  % Convert to microseconds

        % Interpolation: replace NaN values with 0 after interpolation
        R_interp = interp1(t, R, time_delay_us, 'linear', 0);  % [Nz x Nx]

        % Accumulate results
        img_data = img_data + R_interp;

        % --- Plotting ---
        if plot
            imagesc(ypixels*1000, xpixels*1000, abs(img_data));
            axis image;
            xlabel('lateral position (mm)');
            ylabel('depth (mm)');
            title(['Emission ' num2str(ii)]);
            drawnow;
        end
    end
end


%%
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
