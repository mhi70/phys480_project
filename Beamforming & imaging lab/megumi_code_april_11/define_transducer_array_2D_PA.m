% [array] = define_transducer_array_2D([array])
%
% Place a transducer array in the computational grid.
%
% Optional inputs:
%    array : a structure containing transducer array characteristics. If
%    this is not input, default parameters will be used to create this
%    structure:
%       - The default depth is the shallowest possible depth in the
%       grid (just under the space required for the pml).
%       - There is a default element width and spacing based on the
%       emission wavelength. The array is created with as many elements as
%       possible that can be fit into the computational grid.
%
% Global variables (structures) to be defined prior to calling this function:
%   kgrid
%   pulse
%   pml

% Version 1.06
% Authors: Laura Cobus & Manfred Linton
% Last edited: 15-01-25 by Laura Cobus
%
% New for last version:
%
% - trying to get directional sources, to limit the edge effects in PW emission
%       - changed to rect elements instead of line
%       - added definition of array.element.height
%       - ensured that one can make an array of one element spanning the entire grid
%       - changed pml.size*kgrid.dy*2 to pml.size*kgrid.dy
%
% New for this version:
%       - shifts the array down, if array.rot_angles is defined. Otherwise, it will rotate about its center and the ends can be outside of the grid.


function [array2] = define_transducer_array_2D_PA(varargin)

global kgrid
global pulse
global pml

rx = 11; % target position, depth [mm]
ry = 12;

if nargin>0
    array = varargin{1};
else
    array = struct(); % empty structure
end

if isfield(array,'element')
    element = array.element;
else
    element = struct();
end

if ~isfield(element,'bandwidth')
    element.bandwidth = 77;  % [%]
end

% fc = pulse.tone_burst_freq_HF;
% lambda_max = pulse.wavelength_ref*fc/(fc-fc*element.bandwidth/2/100);  % [m]
% lambda_min = pulse.wavelength_ref*fc/(fc+fc*element.bandwidth/2/100);  % [m]
%%% DEFAULT PARAMETERS: Verasonics L11-5v transducer
if ~isfield(element,'pitch') && ~isfield(element,'width')
    %     element.pitch = 0.3e-3; % [m] includes pitch (separation between centres of elements)
    %     element.pitch = round(lambda_max*0.9,1,'significant');
    %     element.pitch = 1e-5*floor(lambda_min*1e5); % or lambda_max, or lambda ?
    element.pitch = 1e-5*floor(pulse.wavelength_ref*1e5)/2;
    element.width = element.pitch;   % [m]
elseif ~isfield(element,'width') && isfield(element,'pitch')
    %     element.width = 0.27e-3;  % [m]
    element.width = element.pitch;  % [m]
else
    % % check that the element isn't larger than the grid
    % if array.element.width>kgrid.y_vec(end)*2 - pml.size*kgrid.dy
    %     disp('   Error! User-defined element size is wider than grid area. ');
    %     disp('   Decreasing to grid area (giving ONE element for the array only).');
    %     disp(' ');
    %     array.element.width=kgrid.y_vec(end)*2-pml.size*kgrid.dy;
    % end
    element.pitch = element.width;   % [m]
end
if ~isfield(element,'height')
    element.height = element.width;
end

y0 = 0;
if isfield(array,'centre')
    y0 = array.centre;
end

if ~isfield(array,'rot_deg')
    array.rot_deg = 0;
end

if ~isfield(element,'num')
    %     element.num = 128;
    array.width=(kgrid.y_vec(end)*2-pml.size*kgrid.dy);
    element.num = floor(array.width/element.pitch);
end
if element.num==0
    element.num=1;
end

pitch_pts = floor(element.pitch/kgrid.dy);

if pml.inside==1
    pml_pts = pml.size*2;
else
    pml_pts=0;
end
if ~isfield(element,'depth')
    depth_in = kgrid.x_vec(1); % [m]
    depth_pts = find(kgrid.x_vec>=depth_in,1) + pml_pts + 1;
    depth_array = kgrid.x_vec(depth_pts);
    element.depth = repmat(depth_array,element.num,1);  % [m]  depth of each element
end

% Example: increase depth of element 50
element.depth(50) = element.depth(50) + rx/1000;  % [m]
if isscalar(element.depth)
    depth_pts = find(kgrid.x_vec>element.depth,1);
    element.depth = kgrid.x_vec(depth_pts); % [m]
else
    if length(element.depth)~=element.num
        disp(' Error! Number of depth points does not match number of elements');
        return
    end
end

% below: putting in user-defined array
array.width = element.pitch*(element.num-1)+element.width;  % [m]

if array.width==element.width  % then there's only one element (pitch is zero)
    element_pos_y_grid = zeros(element.num,1);
else
    % check that the array fits into the grid
    if array.width>kgrid.y_vec(end)*2-pml.size*kgrid.dy
        disp('Error! User-defined array is wider than grid area:');
        disp(['Array: ' num2str(ceil(array.width*1000)) ' mm; Grid area: ' num2str(kgrid.y_vec(end)*2*1000) ' mm.']);
        return
    end
    % below: central positions of each element
    element_pos_y_grid = ((1:element.num) - (element.num+1)/2)*pitch_pts; % [pts]
    element_pos_y_grid = element_pos_y_grid*kgrid.dy + y0; % [m]
end
element.lateral = element_pos_y_grid;

if isfield(array,'rot_angles')==1
    % If array2 is rotated, it gets rotated about x=0.
    % Then, need to set the overall array depth to that of the end element so that the array is still in the grid (otherwise doesn't emit)...

    max_dist = abs(sind(array.rot_angles))*abs(kgrid.y_vec(end)*2)/2;
    element.depth = max(max_dist) + kgrid.x_vec(1);
end

% Create empty array
% karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);   % http://www.k-wave.org/documentation/kWaveArray.php
karray = kWaveArray('BLITolerance', 0.1, 'UpsamplingRate', 10);   % http://www.k-wave.org/documentation/kWaveArray.php

% Add elements to the array
for ind = 1:element.num
    % add element (note: set rotation angle 'theta' to match the global rotation angle)
    % karray.addRectElement([element.depth(ind), element_pos_y_grid(ind)], element.height, element.width, 0); % last argument: theta

    %%% addRectElement(position, Lx, Ly, theta), where
    %     Lx - Height of rect (along x-axis before rotation) [m].
    %     Ly - Width of rect (along y-axis before rotation) [m].
    %       and the x-axis is depth, as defined by kgrid
    % karray.addRectElement([element.depth(ind)-element.height/2, element_pos_y_grid(ind)-element.width/2], element.height, element.width, array.rot_deg); % last argument: theta
    karray.addRectElement([element.depth(ind)-element.height/2, element_pos_y_grid(ind)], element.height, element.width, array.rot_deg); % last argument: theta

    %%% addLineElement(start_point, end_point)
    % karray.addLineElement([element.depth(ind),element_pos_y_grid(ind)-element.width/2], [element.depth(ind),element_pos_y_grid(ind)+element.width/2]); % last argument: theta

    %  http://www.k-wave.org/documentation/kWaveArray.php
    %     addLineElement(start_point, end_point)
    %     start_point	Start coordinate for the line given as a one (1D), two (2D), or three (3D) element vector [m].
    %     end_point	End coordinate for the line given as a one (1D), two (2D), or three (3D) element vector [m].

    % mask = karray.offGridPoints(kgrid, 1, scale)

end

array2.element = element;
array2.karray = karray;

end
