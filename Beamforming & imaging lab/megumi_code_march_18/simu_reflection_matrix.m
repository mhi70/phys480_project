% simu_reflection blank
%
% K-wave code to test acoustic propagation simulations.
%
%%%%%%%%%%%%%%  Notes for running this code:
%kwavekkk
% In script kspaceFirstOrder_initialiseFigureWindow, change img = figure(2); to img = gcf();
%
% In this script, note that x is depth, y is lateral distance

% Version 1.02
%
% Author: Laura Cobus
% Last updated: Sept. 4, 2024
%
% New for this version:
%   - array is defined using default parameters (puts a maximum number of
%   elements into the grid, and at the default depth just under the pml
%   layer). We no longer manually define the elements.


% Edited by Megumi Hirose
% March. 11, 2025
%% ========================================================================
% SET UP, DEFINE PATHS TO FILES & DATA
% =========================================================================

clear variables
clear global
close all hidden
clc

Laura = false;
run_sim = true;

if Laura
    % folder.base = 'C:\Users\lco137\OneDrive - University of Canterbury\Project management\Imaging Team NZ\Megumi 2025\Megumi_photoacoustic_aberration_shared_2025';  % \photoacoustic simulation
    folder.base = 'C:\Users\lco137\OneDrive - University of Canterbury\Project management\Imaging Team NZ\Megumi 2025\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_march_18';
    cd(folder.base)
    addpath(genpath('/Users/travisdunningham/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/k-Wave'));
    addpath([folder.base '/subfunctions']);

    folder.output = '/Users/travisdunningham/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Imaging/Code/Output';
else
    folder.base = 'C:\Users\mhi70\OneDrive - University of Canterbury\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_march_18';
    cd(folder.base)
    folder.kwave = 'C:\Users\mhi70\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\k-Wave';
    addpath(genpath(folder.kwave)); %C:\Users\lcob809\Documents\MATLAB\k-Wave'));
    addpath([folder.base '/subfunctions']); % or wherever the subfunctions folder is
    folder.output ='C:\Users\mhi70\OneDrive - University of Canterbury\Megumi_photoacoustic_aberration_shared_2025\photoacoustic simulation\output';
end
set(0,'DefaultFigureWindowStyle','docked')  % jucst a preference - you can change this if you like


%% ========================================================================
% DEFINE GLOBAL VARIABLES
% =========================================================================

global kgrid
global pml
global medium
global pulse

if run_sim

    %% ========================================================================
    % SET SYSTEM GEOMETRY - 2D
    % =========================================================================
    %   NOTE: For experiments and in mathematical descriptions of beamforming,
    %   we use z to describe position in depth and x as position in the lateral direction.
    %   However, kwave uses x for depth and y for lateral direction, so
    %   some of the variables and definitions coming from this package
    %   (e.g. in the 'kgrid' structure) will use this convention instead.
    %
    %   This can be confusing when coding. It helps to make careful comments (see the lines below).

    yspan = 21; % span in lateral direction [mm]
    xspan = 15; % span depth [mm]

    % =====================================================
    % Create the computational grid. All geometry of these simulations is defined relative to this computational grid.

    % LATERAL
    Ny = 236;     % number of grid points in the y (column) direction - should be odd
    % if rem(Ny,2)==0 % is Ny odd?
    %     Ny = Ny+1;  % if not, make it odd
    % end
    dy = yspan/Ny*1e-3;  % grid point spacing in the y direction [m] (lateral)

    % DEPTH
    Nx =236;
    % if rem(Nx,2)==0
    %     Nx = Nx+1;
    % end
    dx = xspan/Nx*1e-3;  % grid point spacing in the x direction [m] (depth)

    kgrid = kWaveGrid(Nx, dx, Ny, dy);

    % Set the Perfectly Matched Layer properties
    %   see also: http://www.k-wave.org/documentation/example_na_controlling_the_pml.php
    % absorbs waves at the edges
    pml.size  = 10;  % [pts]   default is 20
    pml.alpha = 5;   % default is 2, but that may be too big
    pml.inside = false;


    %%% ========================================================================
    % DEFINE PROPAGATION MEDIUM PROPERTIES
    % =========================================================================

    % homogeneous medium with acoustic wavespeed given by compressional_wavespeed_tissue.
    compressional_wavespeed_tissue = 1540;   % [m/s]
    density_tissue     = 910;  % [kg/m^3]


    % set base medium properties (required to define array in next section)
    medium.density = density_tissue*ones(Nx, Ny);  % [kg/m^3]
    medium.sound_speed = compressional_wavespeed_tissue*ones(Nx, Ny); % [m/s]
    medium.sound_speed_ref = compressional_wavespeed_tissue;

    % define base emission pulse properties
    pulse.tone_burst_freq_HF = 3.25e6;%   % central frequency of input acoustic signal [Hz]
    pulse.tone_burst_cycles_HF = 2; % number of cycles of input pulse
    pulse.wavelength_ref = compressional_wavespeed_tissue/pulse.tone_burst_freq_HF; % wavelength of acoustic pulse travelling through soft tissue [m]
    pulse.magnitude = 1.5e5; % [Pa]

    %%% Check that our geometry and emission pulse are ok for simulation
    mult = 5;
    if(pulse.wavelength_ref<=mult*dy) % check for the lateral direction (y)
        disp(['Error! Wavelength should be >' num2str(mult) ' pts (lateral)']);
        return
    end
    if(pulse.wavelength_ref<=mult*dx) % check for the lateral direction (x)
        disp(['Error! Wavelength should be >' num2str(mult) ' pts (depth)']);
        returna
    end

    % %% Add Noise
    noise = wgn(Nx, Ny, 1); % white gaussian noise added
    noise = (noise/max(noise(:)))*medium.sound_speed_ref*0.15; % scale noise to 15% of sound speed
    medium.sound_speed = medium.sound_speed + noise;


    %% Add a circular reflective target
    rx = 11; % target position, depth [mm]
    % ry = 12; % target position, lateral direction [mm]
    ry = 9;

    r = 0.45*pulse.wavelength_ref*1e3;

    % target radius in multiples of wavelength [mm]
    refl = 3; % reflectivity
    % converting target position in mm to pts on the grid
    rpts = r/dx*1e-3; %[points]
    xc = rx/dx*1e-3; % center of circle x
    yc = ry/dy*1e-3; % center of circle y
    % placing the target into the grid

    for i = 1:Nx
        for j=1:Ny
            dist = sqrt((i-xc).^2+(j-yc).^2);
            if any(dist < rpts)
                medium.sound_speed(i,j) = ...
                    medium.sound_speed(i,j)*refl(dist<rpts);
            end
        end
    end

    %% Plot grid geometry and medium characteristics
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


    %% ========================================================================
    % DEFINE SOURCE TRANSDUCER ARRAY
    % =========================================================================
    % The structure 'array' holds all of the info on the transducer array.
    %   It has a substructure 'element', which holds info in individual array elements

    clear array*

    [array] = define_transducer_array_2D();


    %% Plot transducer array/s on the computational grid
    figure(fig1);
    hold(ax_density,'on')
    plot(ax_density,array.element.lateral*1000,array.element.depth*1000,'wo','MarkerFaceColor','w','MarkerSize',3,'Marker','s');   % plot a line of points where the source array is
    hold(ax_speed,'on')
    plot(ax_speed,array.element.lateral*1000,array.element.depth*1000,'wo','MarkerFaceColor','w','MarkerSize',3,'Marker','s');   % plot a line of points where the source array is


    %% ========================================================================
    % DEFINE ACOUSTIC PULSE FOR EMISSION
    % =========================================================================

    define_input_pulse;

    fig2 = figure;
    plot(pulse.time_us,pulse.signal);
    title('Input pulse - time domain');
    xlabel('$t$ [$\mu$s]','Interpreter','latex');
    ylabel('amplitude');


    %% ========================================================================
    % DEFINE EMISSION
    % =========================================================================

    % % Choose which elements will emit:
    % array.element.emission  = [25];  % [pts]
    %
    % % add time delays
    % array.delays=zeros(1,array.element.num);
    % % array.delays(array.element.emission) =delays;
    %
    % array.nb_sources = length(array.element.emission);
    %
    % disp(' '); disp(['There are ' num2str(array.element.num) ' elements on the emission transducer array']);
    % disp([num2str(array.nb_sources) ' elements will emit.']);
    % disp([num2str(array.element.num) ' elements will record.']);
    % disp(' ');

    %% ========================================================================
    % RUN THE SIMULATION
    % =========================================================================

    % for full matrix capture...

    for kk=1:array.element.num % for every transducer array

        % % Choose which elements will emit
        array.element.emission = kk;

        if kk==1
            [K,array,~] = run_simulation(array,[], true);

            M = zeros(array.element.num, array.element.num, size(K,2));
        else
            [K,array,~] = run_simulation(array,[], false);
        end

        M(kk,:,:) = K;

    end
    %%
    save('Kuu_ver2.mat','M','array','pulse','pml','medium','kgrid','folder'); % save the dataset!


else
    load('Kuu_ver2.mat');
end

K = M;
%% ========================================================================
% Clean up the recorded data
% =========================================================================

K = K - mean(K, 2); % get rid of any constant offset

%% ========================================================================
% VISUALIZE THE RECORDED DATA
% =========================================================================

% Plot the detected signal for various detection array elements.

t = (kgrid.t_array + pulse.length_s/2)*1e6; % time [us]

fig3 = figure;
plot(t, squeeze(K(40,:)));
xlabel('t [us]')
ylabel('signal intensity')
title('signal detected by element 40')

for i=1:array.element.num % loop over all detected signals
    detector_x = array.element.lateral(i); % [m]
    emitter_x = array.element.lateral(array.element.emission(1)); % [m]
    distance = abs(detector_x - emitter_x); % depth is the same for all
    t_direct = pulse.length_s + 2*distance/medium.sound_speed_ref; % travel time straight from emitter to detector [s] could add + pulse period instead instead of multiplying by 2. x 2 doesn't work if detector too close to emitter, e.g 26 where emitter is 25
    npts_direct = round(t_direct*pulse.Fs); % pulse.Fs = samples/second ==> number of samples
    K(i, 1:npts_direct) = 0;
end

fig4 = figure;
scale=500;
hold on;
plot(pulse.time_us, pulse.signal(array.element.emission(1),:), 'DisplayName', 'emitted pulse'); % emitted pulse
plot(t, scale*K(10,:), 'DisplayName', 'detected [10]'); % scaled
hold off;

title('Reflected pulse only')
ylabel('amplitude')
xlabel('t [us]')

legend();

%% Convert K to decibels for plotting
K_dB = 20*log10(abs(K).^2/max(abs(K(:))).^2);

% Make a 2D visualization of all of the recorded signals
x = array.element.lateral*1000; % lateral distance [mm]
fig5 = figure;
imagesc(x, t, abs(squeeze(K_dB)).');
colorbar;
xlabel('lateral distance [mm]');
ylabel('time [us]');


%% ========================================================================
% Look at the frequency content
% =========================================================================
% input pulse
signal = squeeze(pulse.signal(array.element.emission,:));
Nt = pulse.signal_length; % signal length [pts]
Dt = Nt/pulse.Fs*1e6; % signal length [us]
freq = (0:(Nt-1))/Dt; % [MHz]
pulse.signal_fft = abs(fft(signal)); % pulse in freq domain

fig6 = figure;
hold on;
plot(freq,pulse.signal_fft,'ko-','DisplayName', 'input pulse')
xlabel('frequency [MHz]')
ylabel('amplitude')
xlim([0,10]);

% detected signal
signal_detected = squeeze(K(array.element.emission(1),:));
Dt = max(t);
freq = (0:(kgrid.Nt-1))/Dt; % [MHz]
signal_detected_fft = abs(fft(signal_detected));
plot(freq, signal_detected_fft,'DisplayName', 'detected signal (element 1)')

legend();
%% Perform frequency filtering - Hanning window
Img.fc = pulse.tone_burst_freq_HF;
Img.freq_lim = [0.1, 7]; % [MHz]

idxFreq=find(freq>=Img.freq_lim(1) & freq<=Img.freq_lim(2));
Nfkeep = length(idxFreq); % num points to keep

% Hanning Filter on selected frequency
HANN=zeros(kgrid.Nt, 1);
HANN(idxFreq)=hann(Nfkeep);
signal_fft_filt = abs(signal_detected_fft.*HANN);
plot(freq, signal_fft_filt(:,1)*1e18,'DisplayName', 'HANN filtered signal');
plot(freq,HANN*2.2e5,'DisplayName', 'HANN filter');
hold off;
legend();

%% Perform time domain filtering - Tukey window
% Apply Hanning Filter to the entire dataset
Kfft = fft(K,[],2);
Kfft = Kfft.*HANN.';

% Go back into the time domain
Kfilt=ifft(Kfft,[],2);

% Apply time-domain filter
TimeFilt = tukeywin(kgrid.Nt,0.025);

Kfilt_time = Kfilt.*TimeFilt.';


figure(fig4);
hold on;
plot(t, scale*Kfilt(10,:),'g-','DisplayName','filtered [10]')

%% ========================================================================
% Beamforming and imaging
% =========================================================================

% %% A-mode imaging
% % plots of depth-dependent reflected intensity for a particular lateral
% % position x - created by emitting and receiving with the same element
%
% % convert time to depth
% z = t/1e6*medium.sound_speed_ref/2*1000; % depth [mm]
%
% fig7 = figure;
% hold on;
% plot(z, abs(squeeze(K_dB(12,:))));
% % plot(z, abs(K_dB(50,:)));
% % plot(z, abs(K_dB(25,:)));
% % plot(z, abs(K_dB));
% xlabel('$z$ [mm]','Interpreter','latex');
% ylabel('intensity [dB]');
% title('A-mode imaging')




%% B-mode imaging

fig7 = figure;

% target position
z_target = rx - abs(kgrid.x_vec(1)-array.element.depth(1) + array.element.height/2)*1000;
x_target = ry + kgrid.y_vec(1)*1000;

% % plot a circle representing the original pressure source
% source_depth_reltoarray = source_depth -source_radius;  % [mm]
% source_lateral_reltoarray = source_lateral - kgrid.y_vec(end)*1000-source_radius;  % [mm]

plot_circle_x = @(t) r*sin(t) + x_target;
plot_circle_z = @(t) r*cos(t) + z_target;




% Image parameters
zpixels = (-2:0.05:6)/1000; % [m]
xpixels = (-6:0.05:6)/1000; % [m]

Nz = length(zpixels);
Nx = length(xpixels);

img_data = zeros(Nz, Nx);

for ii = 1:array.element.num 

    emitter_x = array.element.lateral(ii); %[m]
    emitter_z = array.element.depth(ii); %[m]

    for xx = 1:array.element.num % for each detection element lateral position
        % location of element
        element_x = array.element.lateral(xx); % [m], element lateral position
        element_z = array.element.depth(xx); %[m], element depth

        R = squeeze(M(ii, xx, :));  % the depth pixels at lateral position xx

        % for each element
        for zz = 1:length(zpixels) % loop over image pixels in depth

            target_z = zpixels(zz); % + kgrid.ky_vec(1)/1000; % [m]

            for xpix = 1:length(xpixels)  % loop over image pixels laterally

                target_x = xpixels(xpix); % [m] same as element lateral position

                transmit_distance = sqrt((emitter_x - target_x)^2+(emitter_z - target_z)^2); % [m], emitter to target
                receive_distance = sqrt((element_x - target_x)^2+(element_z - target_z)^2); % [m], target to element

                time_delay = (transmit_distance + receive_distance)/medium.sound_speed_ref; % [s] total time delay
                time_delay = time_delay * 1e6; % [us]

                % apply time delays to the data
                tmp = interp1(t, R, time_delay); % [us]
                % display(tmp);
                tmp(isnan(tmp))=0; % Avoid NaNs

                img_data(zz,xpix) = img_data(zz,xpix) + tmp;

            end
        end
    % you can watch how the image forms with each subsequent detection:
    imagesc(xpixels*1000, zpixels*1000, abs(img_data))
    axis image
    ylabel('depth (mm)')
    xlabel('lateral position (mm)')
    title(['Emission ' num2str(ii)]);% ', Detection ' num2str(xx)])

    %%% to superimpose target onto image
    hold on;
    fplot(plot_circle_x,plot_circle_z,'k:','LineWidth',1.5);

    drawnow
    end

end




%% B-mode imaging - VECTORIZED

fig7 = figure;

% target position
z_target = rx - abs(kgrid.x_vec(1)-array.element.depth(1) + array.element.height/2)*1000;
x_target = ry + kgrid.y_vec(1)*1000;

% % plot a circle representing the original pressure source
% source_depth_reltoarray = source_depth -source_radius;  % [mm]
% source_lateral_reltoarray = source_lateral - kgrid.y_vec(end)*1000-source_radius;  % [mm]

plot_circle_x = @(t) r*sin(t) + x_target;
plot_circle_z = @(t) r*cos(t) + z_target;




% Image parameters
zpixels = (-2:0.05:6)/1000; % [m]
xpixels = (-6:0.05:6)/1000; % [m]

Nz = length(zpixels);
Nx = length(xpixels);

img_data = zeros(Nz, Nx);

for ii = 1:array.element.num 

    emitter_x = array.element.lateral(ii); %[m]
    emitter_z = array.element.depth(ii); %[m]

    for xx = 1:array.element.num % for each detection element lateral position
        % location of element
        element_x = array.element.lateral(xx); % [m], element lateral position
        element_z = array.element.depth(xx); %[m], element depth

        R = squeeze(M(ii, xx, :));  % the depth pixels at lateral position xx

        % for each element
        for zz = 1:length(zpixels) % loop over image pixels in depth

            target_z = zpixels(zz); % + kgrid.ky_vec(1)/1000; % [m]

            % for xpix = 1:length(xpixels)  % loop over image pixels laterally

                transmit_distance = sqrt((emitter_x - xpixels).^2+(emitter_z - target_z).^2); % [m], emitter to target
                receive_distance = sqrt((element_x - xpixels).^2+(element_z - target_z).^2); % [m], target to element

                time_delay = (transmit_distance + receive_distance)/medium.sound_speed_ref; % [s] total time delay
                time_delay = time_delay * 1e6; % [us]

                % apply time delays to the data
                tmp = interp1(t, R, time_delay); % [us]
                % display(tmp);
                tmp(isnan(tmp))=0; % Avoid NaNs

                img_data(zz,:) = img_data(zz,:) + tmp;

            % end
        end
    % you can watch how the image forms with each subsequent detection:
    imagesc(xpixels*1000, zpixels*1000, abs(img_data))
    axis image
    ylabel('depth (mm)')
    xlabel('lateral position (mm)')
    title(['Emission ' num2str(ii)]);% ', Detection ' num2str(xx)])

    %%% to superimpose target onto image
    hold on;
    fplot(plot_circle_x,plot_circle_z,'k:','LineWidth',1.5);

    drawnow
    end

end


%%
KdB_bmode = 20*log10(abs(img_data).^2/max(abs(img_data(:))).^2);
fig8=figure;
imagesc(x, kgrid.ky_vec, KdB_bmode)
colorbar
clim([-50 0])


