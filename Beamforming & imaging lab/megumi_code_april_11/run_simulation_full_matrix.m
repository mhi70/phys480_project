
% Run simulation for matrix capture simulation
% i.e emit and record from every transducer array
% stores results in 'Kuu.m'

% Edited by Megumi Hirose
% April. 9, 2025


%% set up file path
clear variables
clear global
close all hidden
clc

Laura = true;
run_sim = false; % run if there is no 'Kuu.mat' folder.

if Laura
    % folder.base = 'C:\Users\lco137\OneDrive - University of Canterbury\Project management\Imaging Team NZ\Megumi 2025\Megumi_photoacoustic_aberration_shared_2025';  % \photoacoustic simulation
    folder.base = 'C:\Users\lco137\OneDrive - University of Canterbury\Project management\Imaging Team NZ\Megumi 2025\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_april_11';
    cd(folder.base)
    addpath(genpath('/Users/travisdunningham/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/k-Wave'));
    addpath([folder.base '/../subfunctions']);

    folder.output = '/Users/travisdunningham/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Imaging/Code/Output';
else
    %folder.base='/Users/megumihirose/Desktop/PHYS480_project/AcousticBeamformingLab_scripts';
    folder.base = 'C:\Users\mhi70\OneDrive - University of Canterbury\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_april_11';
    cd(folder.base)
    folder.kwave = 'C:\Users\mhi70\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\k-Wave';
    addpath(genpath(folder.kwave)); %C:\Users\lcob809\Documents\MATLAB\k-Wave'));
    addpath([folder.base '/subfunctions']); % or wherever the subfunctions folder is
    folder.output ='C:\Users\mhi70\OneDrive - University of Canterbury\Megumi_photoacoustic_aberration_shared_2025\photoacoustic simulation\output';
end
set(0,'DefaultFigureWindowStyle','docked')  % jucst a preference - you can change this if you like


%% define global variables

global kgrid
global pml
global medium
global pulse


%% set system geometry (2D)
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
dy = yspan/Ny*1e-3;  % grid point spacing in the y direction [m] (lateral)

% DEPTH
Nx =236;
dx = xspan/Nx*1e-3;  % grid point spacing in the x direction [m] (depth)

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Set the Perfectly Matched Layer properties
pml.size  = 10;  % [pts]   default is 20
pml.alpha = 5;   % default is 2, but that may be too big
pml.inside = false;


%% define propagation medium properties

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
    return
end

%% Add Noise
noise = wgn(Nx, Ny, 1); % white gaussian noise added
noise = (noise/max(noise(:)))*medium.sound_speed_ref*0.15; % scale noise to 15% of sound speed
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


%% define source transducer array
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


%% define acoustic pulse for emission
define_input_pulse;

fig2 = figure;
plot(pulse.time_us,pulse.signal);
title('Input pulse - time domain');
xlabel('$t$ [$\mu$s]','Interpreter','latex');
ylabel('amplitude');

%% run simulation
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
save('Kuu_new.mat','rx','ry','M','array','pulse','pml','medium','kgrid','folder'); % save the dataset!



