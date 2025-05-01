% Created 01 April 2025
% Megumi Hirose
% Rewriting my own Ultra-Sound simulation code 

% reset for each run
clear variables
clear global
close all hidden
clc

folder.base = 'C:\Users\mhi70\OneDrive - University of Canterbury\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_march_18';
% folder.base = 'C:\Users\lco137\OneDrive - University of Canterbury\Project management\Imaging Team NZ\Megumi 2025\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_march_18';
cd(folder.base)
folder.kwave = 'C:\Users\mhi70\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\k-Wave';
addpath(genpath(folder.kwave)); 
addpath('subfunctions'); 
set(0,'DefaultFigureWindowStyle','docked')  % just a preference - you can change this if you like

% Define global variables
global kgrid
global pml
global medium
global pulse

% set system geometry
%   NOTE: For experiments and in mathematical descriptions of beamforming,
%   we use z to describe position in depth and x as position in the lateral direction. 
%   However, kwave uses x for depth and y for lateral direction, so 
%   some of the variables and definitions coming from this package 
%   (e.g. in the 'kgrid' structure) will use this convention instead.

yspan = 21e-3; % span in lateral direction [m]
xspan = 15e-3; % span depth [m]

% Create the computational grid. All geometry of these simulations is defined relative to this computational grid.
% LATERAL
Ny = 256;     % number of grid points in the y (column) direction
dy = yspan/Ny;  % [m] grid point spacing in the y direction [m] (lateral)

% DEPTH
Nx =256;
dx = xspan/Nx;  % [m] grid point spacing in the x direction [m] (depth)

kgrid = kWaveGrid(Nx, dx, Ny, dy); 

% Set the Perfectly Matched Layer properties (absorbs waves at the edges) 
pml.size  = 10;  % [pts]   default is 20
pml.alpha = 5;   % default is 2, but that may be too big
pml.inside = true;


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
pulse.tone_burst_freq_HF = 3.25e6; % [Hz] central frequency of input acoustic signal
pulse.tone_burst_cycles_HF = 2;    % number of cycles of input pulse
pulse.wavelength_ref = compressional_wavespeed_tissue/pulse.tone_burst_freq_HF; % [m] wavelength of acoustic pulse travelling through soft tissue 
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

% Add noise
noise = wgn(Nx, Ny, 1); % white gaussian noise 
noise = (noise/max(noise(:)))*medium.sound_speed_ref*0.15; % scale noise to 15% of sound speed
medium.sound_speed = medium.sound_speed + noise;

% Place reflective target
rx = 11e-3; % target position, depth [m]
ry = 12e-3; % target position, lateral direction [m]

r = 0.45*pulse.wavelength_ref;

% target radius in multiples of wavelength [m]
refl = 3; % reflectivity

% converting target position in m to pts on the grid
rpts = r/dx; %[points]
xc = rx/dx; % center of circle x
yc = ry/dy; % center of circle y

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

clear array* 
[array] = define_transducer_array_2D(); % Define transducer array

define_input_pulse; % Define acoustic pulse for emission

array.element.emission = [1, 10, 50]; % define emitting elements

%% run simulation
[K,array,~] = run_simulation(array,[],true); 

