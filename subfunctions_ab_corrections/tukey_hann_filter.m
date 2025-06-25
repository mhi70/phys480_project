function K_US = tukey_hann_filter(data_file)
load(data_file);
t = (kgrid.t_array + pulse.length_s/2)*1e6; % time [mu s]
Dt = max(t);
freq = (0:(kgrid.Nt-1))/Dt; % [MHz]
for ii = 1:array.element.num % for every emitter
    K = squeeze(K_US(ii,:,:));
    K = K - mean(K, 2); % get rid of any constant offset
    %% Remove emission pulse (only want reflected pulse)
    for kk=1:array.element.num % loop over all detectors
        detector_x = array.element.lateral(kk); % [m]
        emitter_x = array.element.lateral(ii); % [m]
        distance = abs(detector_x - emitter_x); % depth is the same for all
        t_direct = pulse.length_s + 2*distance/medium.sound_speed_ref; % travel time straight from emitter to detector [s]
        npts_direct = round(t_direct*pulse.Fs); % pulse.Fs = samples/second ==> number of samples
        K(kk, 1:npts_direct) = 0;
    end

    %% Perform frequency filtering - Hanning window

    Img.fc = pulse.tone_burst_freq_HF;
    Img.freq_lim = [0.1, 7]; % [MHz]

    idxFreq=find(freq>=Img.freq_lim(1) & freq<=Img.freq_lim(2));
    Nfkeep = length(idxFreq); % num points to keep

    % Hanning Filter on selected frequency
    HANN=zeros(kgrid.Nt, 1);
    HANN(idxFreq)=hann(Nfkeep);

    %% Perform time domain filtering - Tukey window
    % Apply Hanning Filter to the entire dataset
    Kfft = fft(K,[],2);
    Kfft = Kfft.*HANN.';

    % Go back into the time domain
    Kfilt=ifft(Kfft,[],2);

    % Apply time-domain filter
    TimeFilt = tukeywin(kgrid.Nt,0.025);

    Kfilt_time = Kfilt.*TimeFilt.';
    K_US(ii, :, :) = Kfilt_time;

end
end