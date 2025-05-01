% Ultra Sound and Photoacoustic simulation
% Megumi Hirose
% 11 April 2025

clear variables
clear global
close all hidden
clc

%% filepath
Laura = true;

if Laura
    folder.base = 'C:\Users\lco137\OneDrive - University of Canterbury\Project management\Imaging Team NZ\Megumi 2025\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_april_11';
    cd(folder.base)
    addpath(genpath('/Users/travisdunningham/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/k-Wave'));
    addpath([folder.base '/../subfunctions'])
else
    %folder.base='/Users/megumihirose/Desktop/PHYS480_project/AcousticBeamformingLab_scripts';
    folder.base = 'C:\Users\mhi70\OneDrive - University of Canterbury\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_april_11';
    cd(folder.base)
    folder.kwave = 'C:\Users\mhi70\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\k-Wave';
    addpath(genpath(folder.kwave));
    addpath([folder.base '/subfunctions']);
end
set(0,'DefaultFigureWindowStyle','docked')  % jucst a preference - you can change this if you like

%% setup grid
global kgrid
global pml
global medium
global pulse
global source

yspan = 21; % span in lateral direction [mm]
xspan = 15; % span depth [mm]

% computational grid
Ny = 236;     % number of grid points in the y (column) direction
dy = yspan/Ny*1e-3;  % grid point spacing in the y direction [m] (lateral)

Nx = 236;
dx = xspan/Nx*1e-3;  % grid point spacing in the x direction [m] (depth)

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Set the Perfectly Matched Layer properties
%   see also: http://www.k-wave.org/documentation/example_na_controlling_the_pml.php
% absorbs waves at the edges
pml.size  = 10;  % [pts]   default is 20
pml.alpha = 5;   % default is 2, but that may be too big
pml.inside = false;

%% define medium 

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

% Check that our geometry and emission pulse are ok for simulation
mult = 5;
if(pulse.wavelength_ref<=mult*dy) % check for the lateral direction (y)
    disp(['Error! Wavelength should be >' num2str(mult) ' pts (lateral)']);
    return
end
if(pulse.wavelength_ref<=mult*dx) % check for the lateral direction (x)
    disp(['Error! Wavelength should be >' num2str(mult) ' pts (depth)']);
    returna
end

% Add Noise
noise = wgn(Nx, Ny, 1); % white gaussian noise added
noise = (noise/max(noise(:)))*medium.sound_speed_ref*0.15; % scale noise to 15% of sound speed
medium.sound_speed = medium.sound_speed + noise;

medium_US = medium;
medium_PA = medium;

%% Source/reflector settings
source_lateral = 12; % [mm]
source_depth = 11; % [mm]
source_radius = 0.45*pulse.wavelength_ref*1e3; % [mm]
source_pressure = 1.5e5; % [Pa]
source_refl = 3; % [< 3]


%% Add a circular reflective target (for US imaging)
rx = source_depth; % target position, depth [mm]
ry = source_lateral;

r = source_radius;

% target radius in multiples of wavelength [mm]
refl = source_refl; % reflectivity

% converting target position in mm to pts on the grid
rpts = r/dx*1e-3; %[points]
xc = rx/dx*1e-3; % center of circle x
yc = ry/dy*1e-3; % center of circle y

% placing the target into the grid
for i = 1:Nx
    for j=1:Ny
        dist = sqrt((i-xc).^2+(j-yc).^2);
        if any(dist < rpts)
            medium_US.sound_speed(i,j) = ...
                medium_US.sound_speed(i,j)*refl(dist<rpts);
        end
    end
end

% Plot grid geometry and medium characteristics
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

%% define transducer array + pulse

% define source transducer array
% The structure 'array' holds all of the info on the transducer array.
%   It has a substructure 'element', which holds info in individual array elements

clear array*

[array] = define_transducer_array_2D();

% Plot transducer array/s on the computational grid
figure(fig1);
hold(ax_density,'on')
plot(ax_density,array.element.lateral*1000,array.element.depth*1000,'wo','MarkerFaceColor','w','MarkerSize',3,'Marker','s');   % plot a line of points where the source array is
hold(ax_speed,'on')
plot(ax_speed,array.element.lateral*1000,array.element.depth*1000,'wo','MarkerFaceColor','w','MarkerSize',3,'Marker','s');   % plot a line of points where the source array is

define_input_pulse; % define acoustic pulse for emission

%% ULTRASOUND SIMU

medium = medium_US;

%    NEW MEDIUM : add a target
% simulation, imaging

%% PA SIMU

% NEW MEDIUM : add the source at the same place as the target
% simulation, imaging

%% Create source/reflector disc
p_x = round((source_depth*1e-3)/dx);
p_y = round((source_lateral*1e-3)/dy);
p_r = round((source_radius*1e-3)/dx);
disc_1 = makeDisc(kgrid.Nx, kgrid.Ny, p_x, p_y, p_r);

% for full matrix capture...

% [array2] = define_transducer_array_2D_PA();
% array2.element.emission = [50];
% [K, array2,~] = run_simulation(array2, [],true)

medium = medium_PA;

%% Run PA Simulation

% define the source of acoustic waves : here, a circle which expands to create a pressure waves
source.p_mask = disc_1;
source.p0 = source.p_mask*source_pressure;

[K_PA, array, ~] = run_PA_simulation(array, [], true);

%% plot both US and PA images side-by-side
