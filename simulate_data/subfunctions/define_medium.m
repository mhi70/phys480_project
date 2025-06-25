% function define_medium
%
% Defines the base medium (homogeneous) for US/PA simulation
% Default parameters: 
% compressional_wavespeed_tissue = 1540;   % [m/s]
% density_tissue     = 910;  % [kg/m^3]
%
% Megumi Hirose (April 2025)

function define_medium(Nx,dx, Ny,dy)

global medium
global pulse
global kgrid

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