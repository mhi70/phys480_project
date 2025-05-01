% load variables from 'Kuu_filtered.mat' before running

close all hidden
clc

% plot individual + compound images in dB
t = (kgrid.t_array + pulse.length_s/2)*1e6; % time [us]

K = M(1,:,:);
% Convert K to decibels for plotting
K_dB = 20*log10(abs(K).^2/max(abs(K(:))).^2);

% Make a 2D visualization of all of the recorded signals
x = array.element.lateral*1000; % lateral distance [mm]
figure;
imagesc(x, t, abs(squeeze(K_dB)).');
colorbar;
xlabel('lateral distance [mm]');
ylabel('time [us]');
title('emission #1')

% this plots like 168 images, so takes a bit of time
for i = 2:array.element.num
    K_i = M(i, :, :);
    % Convert K to decibels for plotting
    K_dB = 20*log10(abs(K_i).^2/max(abs(K_i(:))).^2);
    
    % Make a 2D visualization of all of the recorded signals
    figure;
    imagesc(x, t, abs(squeeze(K_dB)).');
    colorbar;
    xlabel('lateral distance [mm]');
    ylabel('time [us]');
    title('emission #',i)
    K = K + K_i; % sum up data

    K_dB_sum = 20*log10(abs(K).^2/max(abs(K(:))).^2);
    
    % Make a 2D visualization of all of the recorded signals
    figure;
    imagesc(x, t, abs(squeeze(K_dB_sum)).');
    colorbar;
    xlabel('lateral distance [mm]');
    ylabel('time [us]');
    title('compound',i)
end

%%


%% B-mode imaging

fig7 = figure;

% location of the emitter
emitter = array.element.emission(1);
emitter_x = array.element.lateral(emitter); %[m]
emitter_z = array.element.depth(emitter); %[m]

zpixels = (-2:0.05:6)/1000; % [m]
xpixels = (-6:0.05:6)/1000; % [m]


Nz = length(zpixels);
Nx = length(xpixels);
img_data = zeros(Nz, Nx);

for ii = 1:array.element.num

for xx = 1:array.element.num % for each detection element lateral position
    % location of element
    element_x = array.element.lateral(xx); % [m], element lateral position
    element_z = array.element.depth(xx); %[m], element depth

    R = squeeze(M(ii, xx, :));  % the depth pixels at lateral position xx

    % for each element 
    for zz = 1:length(zpixels) % loop over image pixels in depth
                
        target_z = zpixels(zz); % + kgrid.ky_vec(1)/1000; % [m]

        for xpix = 1:length(xpixels)  % loop over image pixels laterally
    
            target_x = xpixels(xpix); % [m] same as element lateral position

            transmit_distance = sqrt((emitter_x - target_x)^2+(emitter_z - target_z)^2); % [m], emitter to target
            receive_distance = sqrt((element_x - target_x)^2+(element_z - target_z)^2); % [m], target to element
    
            time_delay = (transmit_distance + receive_distance)/medium.sound_speed_ref; % [s] total time delay
            time_delay = time_delay * 1e6; % [us]
    
            % apply time delays to the data
            tmp = interp1(t, R, time_delay); % [us]
            % display(tmp);
            tmp(isnan(tmp))=0; % Avoid NaNs

            img_data(zz,xpix) = img_data(zz,xpix) + tmp;

        end
    end
end
    % you can watch how the image forms with each subsequent detection: 
    imagesc(zpixels*1000, xpixels*1000, abs(img_data))
    axis image
    ylabel('depth (mm)')
    xlabel('lateral position (mm)')
    title(['Emission ' num2str(ii)]);% ', Detection ' num2str(xx)])
    drawnow
end