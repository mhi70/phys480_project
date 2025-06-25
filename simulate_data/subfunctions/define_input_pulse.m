% function define_input_pulse
%
% Defines an input acoustic pulse for kwave numerical simulations
%   - adds and updates some fields in the structures 'pulse' and 'kgrid'
%
% Global variables:
%   kgrid:  structure defining kwave computational grid
%   pulse:  structure defining input pulse 
%   medium: structure defining the medium for the simulation
%
% pulse.Fs is sampling frequency, in units of [pts/second]
% pulse.reflection = 1 (default) indicates a reflection experiment. The detectors
%    record for enough time for the pulse to travel to the farthest point in
%   the computational grid, and back again.  If pulse.reflection = 0, the
%   recorded time is only half of this (use this option for a transmission
%   experiment).

% Version 1.04
% Author: Laura Cobus
% Last updated: 01-dec-23
%
% New for this version:
%   - integrating permanently filterTimeSeries

function define_input_pulse

global kgrid
global pulse
global medium

if isfield(pulse,'reflection')==0
    pulse.reflection = 1;
end

% c is compressional_wavespeed_soft_tissue
c = medium.sound_speed_ref;

% Set the CFL. This controls the stability of the simulation: for more info,
%   go here:   https://www.simscale.com/blog/2017/08/cfl-condition/ 
cfl = 0.15;  % equal or less than 1: the lower the better

tone_burst_freq_HF=pulse.tone_burst_freq_HF;
tone_burst_cycles_HF = pulse.tone_burst_cycles_HF;

if(isfield(pulse,'signal_offset')==0)
    pulse.signal_offset=0;
end

% set end time
t_end = ((pulse.reflection+1))*sqrt((kgrid.x_vec(end)-kgrid.x_vec(1))^2+((kgrid.y_vec(1)-kgrid.y_vec(end)))^2)/c;

% create the time array
kgrid.t_array = makeTime(kgrid, c, cfl, t_end);
pulse.Fs = 1/kgrid.dt;  % temporal sampling frequency [Hz, i.e. samples/sec]

pulse.signal_length = ceil(pulse.Fs/tone_burst_freq_HF*tone_burst_cycles_HF+max(pulse.signal_offset*pulse.Fs)+1); % signal length in time [pts]

pulse.signal = toneBurst(pulse.Fs, tone_burst_freq_HF, tone_burst_cycles_HF,'SignalLength',pulse.signal_length,'SignalOffset',pulse.signal_offset*pulse.Fs);
% filter the input signal - gets rid of high frequencies that can't propagation through the grid (due to spatial sampling restrictions)
for nn=1:size(pulse.signal,1)
    pulse.signal(nn,:) = filterTimeSeries(kgrid, medium, squeeze(pulse.signal(nn,:)),'ZeroPhase',true);%,'PlotSignals',1,'PlotSpectrums',1);
end
pulse.signal = pulse.magnitude*pulse.signal./max(pulse.signal,2);
pulse.time_us = (1:length(pulse.signal))/pulse.Fs*1e6; % [us]

% find the position of pulse peak without delays added - for imaging afterwards
signal_length_min = ceil(pulse.Fs/tone_burst_freq_HF*tone_burst_cycles_HF+1); % signal length in time [pts]
foo = toneBurst(pulse.Fs, tone_burst_freq_HF, tone_burst_cycles_HF,'SignalLength',signal_length_min,'SignalOffset',0);

pulse.length_s =  length(foo)/pulse.Fs;  % [s]

end