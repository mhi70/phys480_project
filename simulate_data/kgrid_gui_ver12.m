% Megumi Hirose (24 April 2025)
% edited 10 June 2025
% New stuff: 
% - photoacoustic simulation is run automatically whenever US simulation is
%   done (both PA and US simulation data is saved to the same file)
% - No noise in the aberration layer
% - seed set to wgn so noise is the same 
% - fixed target so it is perfectly circular (it was an ellipse before)


function kgrid_gui_ver12

% Globals
global kgrid pml medium pulse source

% Create figure
fig = uifigure('Name', 'k-Wave Grid Creator', 'Position', [100, 100, 700, 750]);

%% Panels
% Grid Settings Panel
gridPanel = uipanel(fig, 'Title', 'Grid Settings', 'Position', [20, 570, 660, 150]);
addLabelField(gridPanel, 'x span (mm):', 20, 90, 21, 'xSpanField');
addLabelField(gridPanel, 'Nx:', 340, 80, 512, 'NxField');
addLabelField(gridPanel, 'y span (mm):', 20, 50, 15, 'ySpanField');
addLabelField(gridPanel, 'Ny:', 340, 40, 512, 'NyField');
addLabelField(gridPanel, 'PML size (pts):', 20, 10, 10, 'PMLField');

% Aberration Panel
abPanel = uipanel(fig, 'Title', 'Aberration Layer', 'Position', [20, 440, 660, 110]);
uicheckbox(abPanel, 'Text', 'Add Aberration Layer', 'Position', [20, 50, 200, 22], ...
    'Tag', 'aberrCheck', 'ValueChangedFcn', @(src, event) toggleAberrFields(fig, src.Value));

addLabelField(abPanel, 'Thickness (mm):', 250, 35, 3, 'AberrThickField');
addLabelField(abPanel, 'Coefficient:', 250, 10, 2, 'AberrCoefField');
toggleAberrFields(fig, 'off')

% Target Panel
default_rx = 15; %[mm]
default_ry = 11; %[mm]
default_radius = 0.3; % [mm]
default_reflectivity = 3; % multiple of sound speed in water
targetPanel = uipanel(fig, 'Title', 'Target Properties', 'Position', [20, 300, 660, 120]);
addLabelField(targetPanel, 'Depth rx (mm):', 20, 70, default_rx, 'rxField');
addLabelField(targetPanel, 'Lateral ry (mm):', 340, 70, default_ry, 'ryField');
addLabelField(targetPanel, 'Radius r (mm):', 20, 30, default_radius, 'rField');
addLabelField(targetPanel, 'Reflectivity:', 340, 30, default_reflectivity, 'reflField');

% Buttons
uibutton(fig, 'Text', 'Create kgrid & Plot', 'Position', [30, 120, 200, 30], ...
    'ButtonPushedFcn', @(btn, event) createKgrid(fig));

end


function addLabelField(parent, text, x, y, defaultValue, tag)
uilabel(parent, 'Text', text, 'Position', [x, y, 120, 22]);
uieditfield(parent, 'numeric', 'Position', [x + 130, y, 80, 22], 'Value', defaultValue, 'Tag', tag);
end


function addLabelTextField(parent, text, x, y, defaultValue, tag)
uilabel(parent, 'Text', text, 'Position', [x, y, 120, 22]);
uieditfield(parent, 'text', 'Position', [x + 130, y, 180, 22], 'Value', defaultValue, 'Tag', tag);
end


function toggleAberrFields(fig, state)
tags = {'AberrThickField', 'AberrCoefField'};
for i = 1:length(tags)
    f = fig.findobj('Tag', tags{i});
    f.Visible = state;
end
end


function createKgrid(fig)

global kgrid pml medium pulse source target

% Get values from GUI
pml.size  = getVal(fig, 'PMLField');
pml.alpha = 5;
pml.inside = false;

x_span_mm = getVal(fig, 'xSpanField');
y_span_mm = getVal(fig, 'ySpanField');

if pml.inside
    Nx = getVal(fig, 'NxField');
    Ny = getVal(fig, 'NyField');
else
    Nx = getVal(fig, 'NxField') - 2*pml.size;
    Ny = getVal(fig, 'NyField') - 2*pml.size;
end

x_span = x_span_mm * 1e-3;
y_span = y_span_mm * 1e-3;

dx = x_span / Nx; % [m]
dy = y_span / Ny; % [m]

kgrid = kWaveGrid(Nx, dx, Ny, dy);

define_medium(Nx, dx, Ny, dy);

% add noise
medium_no_noise = medium;
rng(2025); % set seed to make same noise every time
noise = wgn(Nx, Ny, 1);
noise = (noise/max(noise(:))) * medium.sound_speed_ref * 0.15;

medium.sound_speed = medium.sound_speed + noise;

% add aberration layer
if getVal(fig, 'aberrCheck')
    thick_mm = getVal(fig, 'AberrThickField');
    ab_coeff = getVal(fig, 'AberrCoefField');
    thick_idx = round(thick_mm * 1e-3 / dx);

    medium.sound_speed(1:thick_idx, :) = ...
        medium_no_noise.sound_speed(1:thick_idx, :) * ab_coeff;

end

% Reflective target from user inputs
rx = getVal(fig, 'rxField'); % [mm]
ry = getVal(fig, 'ryField'); % [mm]
refl = getVal(fig, 'reflField');
r = getVal(fig, 'rField');   % [mm]

% Save target info
target.dx = rx;  % [mm] (not used here, maybe for display)
target.dy = ry;
target.r = r;

% Convert center and radius to meters
xc = rx * 1e-3;  % [m]
yc = ry * 1e-3;  % [m]
r_m = r * 1e-3;  % [m]

medium_PA = medium;
medium_US = medium;

% Loop through all grid points
for i = 1:Nx
    for j = 1:Ny
        % Physical coordinates of grid point (i, j)
        x = (i - 1) * dx;
        y = (j - 1) * dy;

        % Compute physical distance from circle center
        dist = sqrt((x - xc)^2 + (y - yc)^2);

        % If inside the circle, modify sound speed
        if dist <= r_m
            medium_US.sound_speed(i,j) = medium.sound_speed(i,j) * refl;
        end
    end
end


% Check for existing fig2 and close if valid
if isappdata(fig, 'fig2') && isvalid(getappdata(fig, 'fig2'))
    close(getappdata(fig, 'fig2'));
end

% Create a new figure and store it
fig2 = figure('Name', 'Sound Speed Plot', 'NumberTitle', 'off');
setappdata(fig, 'fig2', fig2);  % Store fig2 handle in fig
imagesc(kgrid.y_vec * 1000, kgrid.x_vec * 1000, medium_US.sound_speed);
title('Sound Speed (m/s)');
xlabel('y: Lateral distance (mm)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('x: Depth (mm)', 'Interpreter', 'latex', 'FontSize', 10);
cb = colorbar;
cb.Label.String = '[m/s]';
axis image;

uibutton(fig, 'Text', 'Run Simulation', 'Position', [250, 120, 200, 30], ...
    'ButtonPushedFcn', @(btn, event) runSimulations(fig, fig2, medium_US, medium_PA));

% Output Panel
outputPanel = uipanel(fig, 'Title', 'Output Settings', 'Position', [20, 200, 660, 80]);
addLabelTextField(outputPanel, 'Output File:', 20, 30, '.mat', 'outputFileField');
end


function val = getVal(fig, tag)
val = fig.findobj('Tag', tag).Value;
end


function runSimulations(fig, fig2, medium_US, medium_PA)
global kgrid pml medium pulse target
% check if file exists
outputFileField = findobj(fig, 'Tag', 'outputFileField');
outputFileName = outputFileField.Value;
if isfile('output/'+""+outputFileName)
    selection = uiconfirm(fig, ['The file "', outputFileName, '" already exists. Overwrite?'], 'File Exists Warning', ...
        'Options', {'Overwrite', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
    if strcmp(selection, 'Cancel')
        return;
    end
end
[array, K_US] = runUS_Simulation(fig, fig2, medium_US);
K_PA = runPA_Simulation(fig, medium_PA);
save(['output/' getVal(fig, 'outputFileField')],'K_US','K_PA','array','medium_US','target','pulse','pml','medium','kgrid'); % save the dataset!
end

function [array, K_US] = runUS_Simulation(fig, fig2,medium_US)
global kgrid pml medium pulse target


medium = medium_US;

% % close set up window
% if ishandle(fig)
%     close(fig);
% end

clear array*

% define source transducer array
[array] = define_transducer_array_2D();

% define acoustic pulse for emission
define_input_pulse;

figure(fig2);
hold on;
plot(array.element.lateral*1000,array.element.depth*1000,'wo','MarkerFaceColor','w','MarkerSize',3,'Marker','s');   % plot a line of points where the source array is

fig3 = figure('Name', 'Input Pulse', 'NumberTitle', 'off');
plot(pulse.time_us,pulse.signal);
title('Input pulse - time domain');
xlabel('$t$ [$\mu$s]','Interpreter','latex');
ylabel('amplitude');

% run simulation
% for full matrix capture...
disp('running')
array.element.emission = 1;

for kk=1:array.element.num % for every transducer array

    % Choose which elements will emit
    array.element.emission = kk;

    if kk==1
        [K,array,~] = run_simulation(array,[], true);

        K_US = zeros(array.element.num, array.element.num, size(K,2));
    else
        [K,array,~] = run_simulation(array,[], false);
    end
    disp([kk,' done'])

    K_US(kk,:,:) = K;

end


end


function K_PA = runPA_Simulation(fig,medium_PA)
global kgrid pml medium pulse target source
medium = medium_PA;
clear array*

% define source transducer array
[array] = define_transducer_array_2D();

% define acoustic pulse for emission
define_input_pulse;

disp('running')

%% Create source/reflector disc
source_pressure = 1.5e5; % [Pa]
p_x = round(target.dx/kgrid.dx*1e-3);
p_y = round(target.dy/kgrid.dy*1e-3);
p_r = round(target.r/kgrid.dx*1e-3);
disc_1 = makeDisc(kgrid.Nx, kgrid.Ny, p_x, p_y, p_r);

%% Run PA Simulation

% define the source of acoustic waves : here, a circle which expands to create a pressure waves
source.p_mask = disc_1;
source.p0 = source.p_mask*source_pressure; % set initial pressure

[K_PA, array, ~] = run_PA_simulation(array, [], true);

end