%% kgrid_gui_ver3 creates user interface to set up kgrid and base medium
% new from ver2: option to add aberrating layer
% Megumi Hirose (16 April 2025)

function kgrid_gui_ver3

global kgrid
global pml
global medium
global pulse

% Create figure
fig = uifigure('Name', '2D k-Wave Grid Creator', 'Position', [100, 100, 620, 580]);

% Default values
default_xspan = 21; % mm
default_yspan = 15; % mm
default_Nx = 256;
default_Ny = 256;
default_PML = 10;
default_aberr_depth = 5;
default_aberr_thickness = 2;
default_aberr_coef = 1.5; % multiple of speed of soundin water

% Inputs
addLabelField(fig, 'x span (mm):', 30, 490, default_xspan, 'xSpanField');
addLabelField(fig, 'Nx:',         320, 490, default_Nx,    'NxField');

addLabelField(fig, 'y span (mm):', 30, 440, default_yspan, 'ySpanField');
addLabelField(fig, 'Ny:',          320, 440, default_Ny,   'NyField');

addLabelField(fig, 'PML size (pts):', 30, 390, default_PML, 'PMLField');

% Aberration layer toggle and options
aberrCheck = uicheckbox(fig, 'Text', 'Add Aberration Layer', 'Position', [30, 340, 200, 22], 'Tag', 'aberrCheck', 'ValueChangedFcn', @(src, event) toggleAberrFields(fig, src.Value));
addLabelField(fig, 'Aberration Depth (mm):', 30, 310, default_aberr_depth, 'AberrDepthField');
addLabelField(fig, 'Aberration Thickness (mm):', 30, 280, default_aberr_thickness, 'AberrThickField');
addLabelField(fig, 'Aberration Coefficient:', 30, 250, default_aberr_coef, 'AberrCoefField');

% Initially hide aberration fields
setAberrationFieldsVisible(fig, false);

% Display dx/dy
uilabel(fig, 'Text', 'dx (mm):', 'Position', [30, 210, 80, 22]);
uilabel(fig, 'Text', 'dy (mm):', 'Position', [320, 210, 80, 22]);
dxLabel = uilabel(fig, 'Text', '', 'Position', [110, 210, 100, 22], 'Tag', 'dxLabel');
dyLabel = uilabel(fig, 'Text', '', 'Position', [400, 210, 100, 22], 'Tag', 'dyLabel');

% Button
uibutton(fig, 'Text', 'Create kgrid & Plot', ...
    'Position', [200, 170, 180, 30], ...
    'ButtonPushedFcn', @(btn, event) createKgrid(fig));

% % Plot area
ax = uiaxes(fig, 'Position', [60, 20, 500, 130]);
xlabel(ax, 'y: lateral (mm)');
ylabel(ax, 'x: depth (mm)');
set(ax, 'Tag', 'PlotAxes');
end

function addLabelField(fig, text, x, y, defaultValue, tag)
uilabel(fig, 'Text', text, 'Position', [x, y, 160, 22]);
uieditfield(fig, 'numeric', 'Position', [x + 170, y, 80, 22], ...
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

    % Get values from GUI
    pml.size  = getVal(fig, 'PMLField');  % [pts]
    pml.alpha = 5;   % default is 2, but that may be too big
    pml.inside = false;

    x_span_mm = getVal(fig, 'xSpanField');
    y_span_mm = getVal(fig, 'ySpanField');
    Nx = getVal(fig, 'NxField') - 2*pml.size;
    Ny = getVal(fig, 'NyField') - 2*pml.size;

    % Convert spans to meters
    x_span = x_span_mm * 1e-3;
    y_span = y_span_mm * 1e-3;

    % Calculate dx, dy
    dx = x_span / Nx;
    dy = y_span / Ny;

    % Display dx/dy in mm
    fig.findobj('Tag', 'dxLabel').Text = num2str(dx * 1e3, '%.4f');
    fig.findobj('Tag', 'dyLabel').Text = num2str(dy * 1e3, '%.4f');

    % Create kgrid
    kgrid = kWaveGrid(Nx, dx, Ny, dy);

    % Define base medium
    define_medium(Nx, dx, Ny, dy);

    % Add aberration if enabled
    if getVal(fig, 'aberrCheck')
        depth_mm = getVal(fig, 'AberrDepthField');
        thick_mm = getVal(fig, 'AberrThickField');
        ab_coeff = getVal(fig, 'AberrCoefField');

        depth_idx = round(depth_mm * 1e-3 / dx);
        thick_idx = round(thick_mm * 1e-3 / dx);

        medium.sound_speed(depth_idx:(depth_idx+thick_idx-1), :) = ...
            medium.sound_speed(depth_idx:(depth_idx+thick_idx-1), :) *ab_coeff;
    end

    % Add Noise
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

    % Plot medium.sound_speed inside GUI
    ax = findall(fig, 'Type', 'axes', 'Tag', 'PlotAxes');
    cla(ax); % clear previous
    imagesc(ax, kgrid.y_vec * 1000, kgrid.x_vec * 1000, medium.sound_speed);
    title(ax, 'Sound Speed (m/s)');
    xlabel(ax, 'y: Lateral distance (mm)', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel(ax, 'x: Depth (mm)', 'Interpreter', 'latex', 'FontSize', 10);
    cb = colorbar(ax);
    cb.Label.String = '[m/s]';
    axis(ax, 'image');
end

function val = getVal(fig, tag)
val = fig.findobj('Tag', tag).Value;
end
