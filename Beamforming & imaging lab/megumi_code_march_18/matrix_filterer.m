% plot matrix M

t = (kgrid.t_array + pulse.length_s/2)*1e6; % time [us]


%% clean recorded data, remove emission pulse
for i = 1:array.element.num
    K_i = squeeze(M(i,:,:));
    K_i = K_i - mean(K_i,2); % remove constant offset
    emitter_x = array.element.lateral(i); % [m]
    for k = 1:array.element.num % for every detector
        detector_x = array.element.lateral(k); % [m]
        distance = abs(detector_x - emitter_x); % depth is the same for all 
        t_direct = pulse.length_s + 2*distance/medium.sound_speed_ref; % travel time straight from emitter to detector [s] could add + pulse period instead instead of multiplying by 2. x 2 doesn't work if detector too close to emitter, e.g 26 where emitter is 25
        npts_direct = round(t_direct*pulse.Fs); % pulse.Fs = samples/second ==> number of samples
        K_i(1:npts_direct) = 0;
    M(i,:,:)=K_i;
    end
end

%%

save('Kuu_cleaned.mat','M','array','pulse','pml','medium','kgrid','folder'); % save the dataset!


%% Perform frequency + time filtering - Hanning + Tukey window

for i = 1:array.element.num
    K_i = squeeze(M(i,:,:)); % data when emitter = i-th element
    signal_detected = squeeze(K_i(i,:));
    Dt = max(t);
    freq = (0:(kgrid.Nt-1))/Dt; % [MHz]
    signal_detected_fft = abs(fft(signal_detected));

    Img.fc = pulse.tone_burst_freq_HF;
    Img.freq_lim = [0.1, 7]; % [MHz]
    
    idxFreq=find(freq>=Img.freq_lim(1) & freq<=Img.freq_lim(2));
    Nfkeep = length(idxFreq); % num points to keep 
    
    % Hanning Filter on selected frequency
    HANN=zeros(kgrid.Nt, 1);
    HANN(idxFreq)=hann(Nfkeep);
    
    % time filter
    TimeFilt = tukeywin(kgrid.Nt,0.025);

    Kfft_i = fft(K_i,[],2); % freq domain
    Kfft_i = Kfft_i.*HANN.'; % apply frequency filter
    Kfilt_i = ifft(Kfft_i,[],2); % time domain
    Kfilt_time = Kfilt_i.*TimeFilt.'; % apply time filter
    M(i,:,:) = Kfilt_time;
end

%%
save('Kuu_filtered.mat','M','array','pulse','pml','medium','kgrid','folder'); % save the dataset!


