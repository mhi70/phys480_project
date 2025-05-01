%% kgrid_gui_ver9 creates user interface to set up kgrid and base medium
% New from ver 8: 
% can change radius of target, 'r'
% saves simulation data
% added PA imaging
%
% Megumi Hirose (24 April 2025)

function kgrid_gui_ver9

global kgrid
global pml
global medium
global pulse
global target

% Create figure
fig = uifigure('Name', 'k-Wave Grid Creator', 'Position', [100, 100, 620, 720]);

% Default values
default_xspan = 21; % [mm]
default_yspan = 15; % [mm]
default_Nx = 256;
default_Ny = 256;
default_PML = 10;
default_aberr_depth = 5; % [mm]
default_aberr_thickness = 2; % [mm]
default_aberr_coef = 1.5; % multiple of speed of sound in water
default_rx = 11; % [mm]
default_ry = 12; % [mm]
default_r =  0.45 * 4.7385e-04 * 1e3; % [mm]
default_refl = 3; % <= 3


% Inputs
addLabelField(fig, 'x span (mm):', 30, 560, default_xspan, 'xSpanField');
addLabelField(fig, 'Nx:',         320, 560, default_Nx,    'NxField');

addLabelField(fig, 'y span (mm):', 30, 510, default_yspan, 'ySpanField');
addLabelField(fig, 'Ny:',          320, 510, default_Ny,   'NyField');

addLabelField(fig, 'PML size (pts):', 30, 460, default_PML, 'PMLField');

% Aberration layer toggle and options
aberrCheck = uicheckbox(fig, 'Text', 'Add Aberration Layer', 'Position', [30, 420, 200, 22], 'Tag', 'aberrCheck', 'ValueChangedFcn', @(src, event) toggleAberrFields(fig, src.Value));
addLabelField(fig, 'Aberration Depth (mm):', 30, 390, default_aberr_depth, 'AberrDepthField');
addLabelField(fig, 'Aberration Thickness (mm):', 30, 360, default_aberr_thickness, 'AberrThickField');
addLabelField(fig, 'Aberration Coefficient:', 30, 330, default_aberr_coef, 'AberrCoefField');

% Initially hide aberration fields
setAberrationFieldsVisible(fig, false);

% Reflective target inputs
addLabelField(fig, 'Target depth rx (mm):', 30, 290, default_rx, 'rxField');
addLabelField(fig, 'Target lateral ry (mm):', 30, 260, default_ry, 'ryField');
addLabelField(fig, 'Target radius r (mm):', 30, 230, default_r, 'rField');
addLabelField(fig, 'Target reflectivity:', 30, 200, default_refl, 'reflField');

% Display dx/dy
uilabel(fig, 'Text', 'dx (mm):', 'Position', [30, 190, 80, 22]);
uilabel(fig, 'Text', 'dy (mm):', 'Position', [320, 190, 80, 22]);
dxLabel = uilabel(fig, 'Text', '', 'Position', [110, 190, 100, 22], 'Tag', 'dxLabel');
dyLabel = uilabel(fig, 'Text', '', 'Position', [400, 190, 100, 22], 'Tag', 'dyLabel');

% Buttons
uibutton(fig, 'Text', 'Create kgrid & Plot', ...
    'Position', [200, 150, 180, 30], ...
    'ButtonPushedFcn', @(btn, event) createKgrid(fig));

end

function addLabelField(fig, text, x, y, defaultValue, tag)
uilabel(fig, 'Text', text, 'Position', [x, y, 160, 22]);
uieditfield(fig, 'numeric', 'Position', [x + 170, y, 80, 22], ...
    'Value', defaultValue, 'Tag', tag);
end


function addLabelTextField(fig, text, x, y, defaultValue, tag)
    uilabel(fig, 'Text', text, 'Position', [x, y, 160, 22]);
    uieditfield(fig, 'text', 'Position', [x + 170, y, 180, 22], ...
        'Value', defaultValue, 'Tag', tag);
end


function setAberrationFieldsVisible(fig, visible)
visibility = 'off';
if visible
    visibility = 'on';
end
fig.findobj('Tag', 'AberrDepthField').Visible = visibility;
fig.findobj('Tag', 'AberrThickField').Visible = visibility;
fig.findobj('Tag', 'AberrCoefField').Visible = visibility;
fig.findobj('Text', 'Aberration Depth (mm):').Visible = visible;
fig.findobj('Text', 'Aberration Thickness (mm):').Visible = visible;
fig.findobj('Text', 'Aberration Coefficient:').Visible = visible;
end

function toggleAberrFields(fig, state)
setAberrationFieldsVisible(fig, state);
end

function createKgrid(fig)

global kgrid
global pml
global medium
global pulse
global target

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

fig.findobj('Tag', 'dxLabel').Text = num2str(dx * 1e3, '%.4f');
fig.findobj('Tag', 'dyLabel').Text = num2str(dy * 1e3, '%.4f');

kgrid = kWaveGrid(Nx, dx, Ny, dy);

define_medium(Nx, dx, Ny, dy);

if getVal(fig, 'aberrCheck')
    depth_mm = getVal(fig, 'AberrDepthField');
    thick_mm = getVal(fig, 'AberrThickField');
    ab_coeff = getVal(fig, 'AberrCoefField');

    depth_idx = round(depth_mm * 1e-3 / dx);
    thick_idx = round(thick_mm * 1e-3 / dx);

    medium.sound_speed(depth_idx:(depth_idx+thick_idx-1), :) = ...
        medium.sound_speed(depth_idx:(depth_idx+thick_idx-1), :) * ab_coeff;
end

noise = wgn(Nx, Ny, 1);
noise = (noise/max(noise(:))) * medium.sound_speed_ref * 0.15;
medium.sound_speed = medium.sound_speed + noise;

% Reflective target from user inputs
rx = getVal(fig, 'rxField'); % [mm]
ry = getVal(fig, 'ryField'); % [mm]
refl = getVal(fig, 'reflField');
r = getVal(fig, 'rField'); % [mm]

target.rx = rx;% [mm]
target.ry = ry;% [mm]
target.refl = refl;
target.r = r;% [mm]

rpts = r/dx*1e-3;
xc = rx/dx*1e-3;
yc = ry/dy*1e-3;

medium_PA = medium;
medium_US = medium;

for i = 1:Nx
    for j = 1:Ny
        dist = sqrt((i - xc)^2 + (j - yc)^2);
        if dist < rpts
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

uibutton(fig, 'Text', 'Run US Simulation', ...
    'Position', [200, 100, 180, 30], ...
    'ButtonPushedFcn', @(btn, event) runUS_Simulation(fig,fig2,medium_US));


uibutton(fig, 'Text', 'Run PA Simulation', ...
    'Position', [400, 100, 180, 30], ...
    'ButtonPushedFcn', @(btn, event) runPA_Simulation(fig,fig2, medium_PA));


addLabelTextField(fig, 'US Output File:', 30, 160, 'K_US.mat', 'usFileField');
addLabelTextField(fig, 'PA Output File:', 30, 130, 'K_PA.mat', 'paFileField');
end


function val = getVal(fig, tag)
val = fig.findobj('Tag', tag).Value;
end


function runUS_Simulation(fig, fig2,medium_US)
global kgrid
global pml
global medium
global pulse

medium = medium_US;
% close set up window
if ishandle(fig)
    close(fig);
end

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

        M = zeros(array.element.num, array.element.num, size(K,2));
    else
        [K,array,~] = run_simulation(array,[], false);
    end
    disp(kk,' done')

    M(kk,:,:) = K;

end
save('output/K_US.mat','M','array','pulse','pml','medium','kgrid'); % save the dataset!

end


function runPA_Simulation(fig, fig2, medium_PA)
global kgrid
global pml
global medium
global pulse
global target
global source
medium = medium_PA;
% close set up window
if ishandle(fig)
    close(fig);
end

clear array*

% define source transducer array
[array] = define_transducer_array_2D();

% define acoustic pulse for emission
define_input_pulse;

disp('running')

%% Create source/reflector disc
source_pressure = 1.5e5; % [Pa]
p_x = round(target.rx/kgrid.dx*1e-3);
p_y = round(target.ry/kgrid.dy*1e-3);
p_r = round(target.r/kgrid.dx*1e-3);
disc_1 = makeDisc(kgrid.Nx, kgrid.Ny, p_x, p_y, p_r);

%% Run PA Simulation

% define the source of acoustic waves : here, a circle which expands to create a pressure waves
source.p_mask = disc_1;
source.p0 = source.p_mask*source_pressure; % set initial pressure

[K_PA, array, ~] = run_PA_simulation(array, [], true);
save('output/K_PA.mat','K_PA','array','medium_PA','target','pulse','pml','medium','kgrid')
end

