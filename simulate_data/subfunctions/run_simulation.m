% function [K, array, array2] = run_simulation(array_source, [array_sensor, PlotSim, ascale])
%
% This function runs the kwave simulation.
%
% Inputs:
%   - array_source
%      You need to define:
%        array_source.element.emission : which array elements will emit (otherwise none will emit)
%        optional:  array_source.delays : what time delays to apply to the emitted signals - the default is 0 [s]
%
% Optional inputs:
%   - array_sensor : optional second array to record (if not defined, the
%   source array will record)
%   - PlotSim : = true for a visualization of the simulation (default).
%               = false to not show the simulation (faster)
%   - ascale :  Controls the amplitude scaling of the simulation.
%               = 'constant' for a manually-defined scaling (default). This is best for basic visualization, as it avoids
%                 the non-physical oscillations at late times when the pulse has left the computational grid.
%               = 'auto' for autoscaling at each time. This may be required when using more complex inputs
%                 (like plane-wave simulations)
%
% Ouputs:
%   - K(emission element, detection element, time) holds the recorded signals

% Version 1.05
% Author: Laura Cobus
% Last updated: 20-jan-2025
%
%   Each sensor records the time-dependent pressure at the sensor point.
%   These recorded signals are held in the variable 'sensor_data'.
%   For each array element, we average the recorded signal for all points in the computational grid that make up that array element
%
% New for this version:
% - can rotate the emission array

function [K, array, array2] = run_simulation(array, varargin)

global kgrid
global medium
global pulse
global pml

disp('Running simulation!...');

% defaults for the pml:
if ~isfield(pml,'size')
    pml.size = 40;
end
if ~isfield(pml,'alpha')
    pml.alpha = 1;
end
if ~isfield(pml,'inside')
    pml.inside = false;
end
% defaults:
array2 = array;
PlotSim = true;
ascale = 'constant';
% optional inputs
if nargin>1 && ~isempty(varargin{1})
    array2 = varargin{1};
end
if nargin>2 && ~isempty(varargin{2})
    PlotSim = varargin{2};
end
if nargin>3 && ~isempty(varargin{3})
    ascale = varargin{3};
end

Nelem = array.element.num;
Nelem2 = array2.element.num;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(array,'rot_deg')
%%% optional:  rotating the emission array

array.karray.setArrayPosition([0,0], array.rot_deg); % rotation is counter-clockwise
%     setArrayPosition(translation, rotation)
%
%         Sets the property array_transformation (an affine transform)
%         based on the values for translation and rotation. The
%         translations are given as [dx, dy] in 2D and [dx, dy, dz] in 3D.
%         The rotations angle/s are given as [th] in 2D (counter-clockwise)
%         and [x_th, y_th, z_th] in 3D (rotation about x then y' then z'').
% 
%         translation - Array translation [m].
%         rotation    - Array rotation [degrees].

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isfield(array,'delays')
    delays = array.delays;
else
    delays = zeros(1,Nelem);
end
if ~isfield(array,'nb_sources')
    array.nb_sources = length(array.element.emission);
end

if ~isfield(array,'up')
    array.up = 'p'; % default
end
if ~isfield(array2,'up')
    array2.up = array.up; % default: same as emission
end

sensor.record = {array2.up};

% disp('Calculating grid masks...');
if array.up=='u'
    source.u_mask = array.karray.getArrayBinaryMask(kgrid);
else
    source.p0 = [];
    source.p_mask = array.karray.getArrayBinaryMask(kgrid);
end

sensor.mask = array2.karray.getArrayBinaryMask(kgrid);
% sensor.frequency_response = [pulse.tone_burst_freq_HF, array.element.bandwidth]; % [fc (Hz), bandwidth (%)]; % only works if you're measuring pressure (not velocity)

% array.element.directivity = pi;
% sensor = add_element_directivity(sensor,array);

% % define the directionality of the sensor elements
% sensor.directivity_angle = zeros(kgrid.Nx, kgrid.Ny);
% sensor.directivity_angle(1, :) = 0;    	 % max sensitivity in x direction
% sensor.directivity_angle(end, :) = 0;  	 % max sensitivity in x direction
% sensor.directivity_angle(:, 1) = pi/2;   % max sensitivity in y direction
% sensor.directivity_angle(:, end) = pi/2; % max sensitivity in y direction
% 
% % define the directivity size
% sensor.directivity_size = 20 * kgrid.dx;

K = zeros(Nelem2,kgrid.Nt,Nelem);

% apply time delays to output pulse
pulse.signal_offset = abs(delays); % [s] -> required for define_input_pulse.m
define_input_pulse;

% only emit from desired elements
emask = zeros(size(pulse.signal));
emask(array.element.emission,:) = 1;
pulse.signal = pulse.signal.*emask;

if array.up=='u'
    source.ux = pulse.signal;
    source.ux = array.karray.getDistributedSourceSignal(kgrid, source.ux);
else
    source.p = pulse.signal;  % [elements, time]    
    source.p = array.karray.getDistributedSourceSignal(kgrid, source.p);
    if isfield(array,'apod_transmit')==1
        if strcmp(array.apod_transmit,'tukey')==1
            source.p = tukeywin(size(source.p,1),0.5).*source.p;
        end
    end    
end

% some options for the simulation and the simulation animation
%     max_pressure_display = max(pulse.signal)/2;%.^2*sqrt(array.nb_sources);
% max_pressure_display = max(source.p(:))/10/sqrt(array.nb_sources)*(array.element.height/array.element.width+1);
% the /10 is to have a pleasing-looking scale. the /sqrt(array.nb_sources) is needed when emitting w/ >1 source
max_pressure_display = 0.5*max(source.p(:))/sqrt(array.nb_sources)*(array.element.height/array.element.width+1);
pdisp = [-max_pressure_display,max_pressure_display]; % for input to simulation
if strcmp(ascale,'constant')
    ascale = pdisp;
end

input_args = {'PMLSize', pml.size, 'PlotSim', PlotSim,'PlotScale',ascale,'PlotPML',true,'PMLAlpha',pml.alpha,'PMLInside',pml.inside, 'DataCast','single'};  % or set PlotScale to 'auto'

% running the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}); %, 'PlotLayout', true, 'PlotPML', false);

% retrieving the data from the simulation
if array2.up=='u'
    K = array2.karray.combineSensorData(kgrid, sensor_data.ux);
else
    K = array2.karray.combineSensorData(kgrid, sensor_data.p);
end

end

