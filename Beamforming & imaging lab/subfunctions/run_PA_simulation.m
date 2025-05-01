% function [K, array, array2] = run_simulation(array_source, [array_sensor, PlotSim])
%
% This function runs the kwave simulation.
%
% You need to define:
%   array_source.element.emission : which array elements will emit
%   array_source.delays : what time delays to apply to the emitted signals - the default is 0 [s]
%
%   Each sensor records the time-dependent pressure at the sensor point.
%   These recorded signals are held in the variable 'sensor_data'.
%   For each array element, we average the recorded signal for all points in the computational grid that make up that array element
%
% K holds the recorded signals.
% It has dimensions:   K(emission element, detection element, time)

% Version 1.01
% Author: Travis Dunningham & Laura Cobus
% Last updated: 26 march 2025 by Laura Cobus
%
% Changelog
% - disabled automatic smoothing of source.p0 in kwave simu (changes shape of array). This is now done manually in this function just before.

function [K, array, array2] = run_PA_simulation(array, varargin)  %medium, pulse, array, PlotSim, varargin)

global kgrid
global medium
global pml
global source

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
if nargin>1
    if ~isempty(varargin{1})
        array2 = varargin{1};
    else
        array2 = array;
    end
    PlotSim = varargin{2};
end

Nelem = array.element.num;
Nelem2 = array2.element.num;



if ~isfield(array,'up')
    array.up = 'p'; % default
end
if ~isfield(array2,'up')
    %array2.up = array.up; % default: same as emission
    array2.up = 'p'; % default: same as emission
end

sensor.record = {array2.up};

%{
% disp('Calculating grid masks...');
if array.up=='u'
    source.u_mask = array.karray.getArrayBinaryMask(kgrid);
else
    source.p0 = [];
    source.p_mask = array.karray.getArrayBinaryMask(kgrid);
end
%}

%% Defines source mask for photoacoustic source (TD)
% disp('Creating source mask');
%if array.up=='u'
%    source.u_mask = zeros(Nx,Ny);
%    source.u_mask(array.p_x,array.p_y) = 1;
%else
%    source.p0 = [];
%    source.p_mask = zeros(kgrid.Nx,kgrid.Ny);
%    source.p_mask(array.p_x,array.p_y) = 1;
%end

% source.p0 = smooth(source.p0);

%%
sensor.mask = array2.karray.getArrayBinaryMask(kgrid);
% sensor.frequency_response = [pulse.tone_burst_freq_HF, array.element.bandwidth]; % [fc (Hz), bandwidth (%)]; % only works if you're measuring pressure (not velocity)

% sensor = add_element_directivity(sensor,array2);

% K = zeros(Nelem2,kgrid.Nt,Nelem);

% apply time delays to output pulse
%pulse.signal_offset = abs(delays); % [s]
%define_PA_pulse;

% only emit from desired elements
%emask = zeros(size(pulse.signal));
%emask(array.element.emission,:) = 1;
%pulse.signal = pulse.signal.*emask; 

%if array.up=='u'
    %source.ux = pulse.signal;
    %source.ux = array.karray.getDistributedSourceSignal(kgrid, source.ux);
%else
   % source.p = pulse.signal;
   % source.p = array.karray.getDistributedSourceSignal(kgrid, source.p);
%end

% some options for the simulation and the simulation animation
%     max_pressure_display = max(pulse.signal)/2;%.^2*sqrt(array.nb_sources);

max_pressure_display = max(source.p0(:))/2;%.^2;%*sqrt(array.nb_sources);
pdisp = [-max_pressure_display,max_pressure_display]; % for input to simulation

input_args = {'PMLSize', pml.size, 'PlotSim', PlotSim,'PlotScale','auto','PlotPML',true,'PMLAlpha',pml.alpha,'PMLInside',pml.inside, 'DataCast','single','Smooth',false};  % or set PlotScale to 'auto'

% running the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}); %, 'PlotLayout', true, 'PlotPML', false);

% retrieving the data from the simulation
if array2.up=='u'
    K = array2.karray.combineSensorData(kgrid, sensor_data.ux);
else
    K = array2.karray.combineSensorData(kgrid, sensor_data.p);
end

end

