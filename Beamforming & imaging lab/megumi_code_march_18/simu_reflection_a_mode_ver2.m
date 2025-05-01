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

% close all hidden
% clc

% below: the folder containing this script
% folder.base = 'P:\2025\MAPH480\AcousticBeamformingLab_scripts';
folder.base = 'C:\Users\lco137\OneDrive - University of Canterbury\Project management\Imaging Team NZ\Megumi 2025\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_march_18';

cd(folder.base)

folder.kwave = 'C:\Users\mhi70\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\k-Wave';
addpath(genpath(folder.kwave)); %C:\Users\lcob809\Documents\MATLAB\k-Wave'));
addpath('subfunctions'); % or wherever the subfunctions folder is

set(0,'DefaultFigureWindowStyle','docked')  % just a preference - you can change this if you like

%% ========================================================================
% DEFINE GLOBAL VARIABLES
% =========================================================================

global kgrid
global pml
global medium
global pulse


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
Ny = 225;     % number of grid points in the y (column) direction - should be odd
if rem(Ny,2)==0 % is Ny odd?
    Ny = Ny+1;  % if not, make it odd
end
dy = yspan/Ny*1e-3;  % grid point spacing in the y direction [m] (lateral)

% DEPTH
Nx =225;
if rem(Nx,2)==0
    Nx = Nx+1;
end
dx = xspan/Nx*1e-3;  % grid point spacing in the x direction [m] (depth)

kgrid = kWaveGrid(Nx, dx, Ny, dy);  % number of grid points in the x (row) direction - should be odd

% Set the Perfectly Matched Layer properties
%   see also: http://www.k-wave.org/documentation/example_na_controlling_the_pml.php
% absorbs waves at the edges 
pml.size  = 10;  % [pts]   default is 20
pml.alpha = 2;   % default is 2, but that may be too big
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
    return
end

%% Add Noise
noise = wgn(Nx, Ny, 1); % white gaussian noise added
noise = (noise/max(noise(:)))*medium.sound_speed_ref*0.05; % scale noise to 15% of sound speed
medium.sound_speed = medium.sound_speed + noise;


%% Add a circular reflective target
rx = 11; % target position, depth [mm]
ry = 12; % target position, lateral direction [mm]

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
plot(pulse.time_us,pulse.signal,'o');
title('Input pulse - time domain');
xlabel('$t$ [$\mu$s]','Interpreter','latex');
ylabel('amplitude');


%% ========================================================================
% DEFINE EMISSION
% =========================================================================

% Choose which elements will emit:

xo = ry*1e-3; % target lateral [m]
xo = xo + kgrid.y_vec(1); % shift  [m]
emitter = find(array.element.lateral*1000>=xo,1);

array.element.emission  = [emitter];  % [pts]

array.nb_sources = length(array.element.emission);

disp(' '); disp(['There are ' num2str(array.element.num) ' elements on the emission transducer array']);
disp([num2str(array.nb_sources) ' elements will emit.']);
disp([num2str(array.element.num) ' elements will record.']);
disp(' '); 

%% ========================================================================
% RUN THE SIMULATION
% =========================================================================
[K,array,~] = run_simulation(array);

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
    t_direct = 2*distance/medium.sound_speed_ref; % travel time straight from emitter to detector [s] could add + pulse period instead instead of multiplying by 2. x 2 doesn't work if detector too close to emitter, e.g 26 where emitter is 25
    npts_direct = round(t_direct*pulse.Fs); % pulse.Fs = samples/second ==> number of samples
    K(i, 1:npts_direct) = 0;
end

fig4 = figure;
scale=500;
hold on;
plot(pulse.time_us, pulse.signal(array.element.emission(1),:), 'DisplayName', 'emitted pulse'); % emitted pulse
% plot(t, scale*K(5,:), 'DisplayName', 'detected [5]'); % scaled
plot(t, scale*K(10,:), 'DisplayName', 'detected [10]'); % scaled
% plot(t, 500*K(25,:), 'DisplayName', 'detected [25]');
% plot(t, K(26,:), 'DisplayName', 'detected [26]');
% plot(t, K(30,:), 'DisplayName', 'detected [30]');
% plot(t, scale*K(50,:), 'DisplayName', 'detected [50]');
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


% add calculated time delay
xo = ry*1e-3; % target lateral [m]
z = rx*1e-3; % target depth [m]

xd = array.element.lateral; % detector lateral  [m]
xe = array.element.lateral(25); % emitter lateral [m]
xo = xo + kgrid.y_vec(1); % shift  [m]

yd = abs(kgrid.x_vec(1)-array.element.depth(1));  % [m]
ye = abs(kgrid.x_vec(1)-array.element.depth(25));  % [m]

v = medium.sound_speed_ref;  % [m/s]

d1 = sqrt((z-ye)^2+(xo-xe)^2); % distance, object to emitter [m]
d2 = sqrt((z-yd)'.^2+(xo-xd).^2); % distance, object to detector [m]
t_delay = (d1+d2)/v*1e6; % [us]


% hold on;
% plot(x, t_delay,'r-')

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
Img.freq_lim = [1, 5]; % [MHz]

Kfilt = FreqFilt_Hanning(K,t,Img.freq_lim,1);

figure(fig4);
hold on;
plot(t, scale*Kfilt(emitter,:),'g-','DisplayName','filtered [10]')

%% ========================================================================
% Beamforming and imaging
% =========================================================================

K_dB = 20*log10(abs(Kfilt).^2/max(abs(Kfilt(:))).^2);

%% A-mode imaging
% plots of depth-dependent reflected intensity for a particular lateral
% position x - created by emitting and receiving with the same element

% convert time to depth 
z = t/1e6*medium.sound_speed_ref/2*1000; % depth [mm]

fig6 = figure;
hold on;
plot(z, (squeeze(K_dB(emitter,:))));
% plot(z, abs(K_dB(50,:)));
% plot(z, abs(K_dB(25,:)));
% plot(z, abs(K_dB));
xlabel('$z$ [mm]','Interpreter','latex');
ylabel('intensity [dB]');
title('A-mode imaging')



