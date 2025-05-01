
function kgrid_gui
global kgrid
global pml
global medium
global pulse

% Create figure
fig = uifigure('Name', '2D k-Wave Grid Creator', 'Position', [100, 100, 620, 520]);

% Default values
default_xspan = 21; % mm
default_yspan = 15; % mm
default_Nx = 256;
default_Ny = 256;
default_PML = 10;

% Inputs
addLabelField(fig, 'x span (mm):', 30, 440, default_xspan, 'xSpanField');
addLabelField(fig, 'Nx:',         240, 440, default_Nx,    'NxField');

addLabelField(fig, 'y span (mm):', 30, 390, default_yspan, 'ySpanField');
addLabelField(fig, 'Ny:',          240, 390, default_Ny,   'NyField');

addLabelField(fig, 'PML size (pts):', 30, 340, default_PML, 'PMLField');

% Display dx/dy
uilabel(fig, 'Text', 'dx (mm):', 'Position', [30, 290, 80, 22]);
uilabel(fig, 'Text', 'dy (mm):', 'Position', [240, 290, 80, 22]);
dxLabel = uilabel(fig, 'Text', '', 'Position', [110, 290, 100, 22], 'Tag', 'dxLabel');
dyLabel = uilabel(fig, 'Text', '', 'Position', [320, 290, 100, 22], 'Tag', 'dyLabel');

% Button
uibutton(fig, 'Text', 'Create kgrid & Plot', ...
    'Position', [200, 250, 180, 30], ...
    'ButtonPushedFcn', @(btn, event) createKgrid(fig));

% Plot area
ax = uiaxes(fig, 'Position', [60, 20, 500, 210]);
title(ax, '2D Grid with PML');
xlabel(ax, 'X (mm)');
ylabel(ax, 'Y (mm)');
set(ax, 'Tag', 'PlotAxes');
end

function addLabelField(fig, text, x, y, defaultValue, tag)
uilabel(fig, 'Text', text, 'Position', [x, y, 120, 22]);
uieditfield(fig, 'numeric', 'Position', [x + 110, y, 80, 22], ...
    'Value', defaultValue, 'Tag', tag);
end

function createKgrid(fig)
global kgrid
global pml
global medium
global pulse
% Get values from GUI

% PML layer changes size of grid, adjust so Nx and Ny are 2^n
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

% define base medium
define_medium(Nx, dx, Ny, dy)
fig1 = figure;
imagesc(kgrid.y_vec*1000,kgrid.x_vec*1000,medium.sound_speed); % plot a 2D image of the medium sound speed
title('sound speed')
xlabel('Lateral distance (mm)', 'Interpreter','latex','Fontsize',10)
ylabel('Depth (mm)', 'Interpreter','latex','Fontsize',10)
a = colorbar;
a.Label.String = '[m/s]'; % units for colorbar
axis image

end

function val = getVal(fig, tag)
val = fig.findobj('Tag', tag).Value;
end





% disp('Created 2D kgrid:');
% disp(kgrid);
% 
% % Plot
% ax = findobj(fig, 'Type', 'uiaxes', 'Tag', 'PlotAxes');
% cla(ax);
% hold on;
% 
% % Full domain
% x_mm = (0:Nx-1) * dx * 1e3;
% y_mm = (0:Ny-1) * dy * 1e3;
% 
% % Create blank grid image
% grid_data = zeros(Ny, Nx);
% 
% % Display grid
% imagesc(ax, x_mm, y_mm, grid_data);
% colorbar(ax, 'Ticks', []);
% axis(ax, 'image');
% xlim(ax, [0 x_span_mm]);
% ylim(ax, [0 y_span_mm]);
% 
% title(ax, 'kgrid');
% hold off;