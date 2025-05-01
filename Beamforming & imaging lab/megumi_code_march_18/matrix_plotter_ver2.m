% load variables from 'Kuu_filtered.mat' before running

close all hidden
clc

% plot individual + compound images in dB
t = (kgrid.t_array + pulse.length_s/2)*1e6; % time [us]

K = M(1,:,:);
% Convert K to decibels for plotting
K_dB = 20*log10(abs(K).^2/max(abs(K(:))).^2);
M(1, :,:) = K_dB;

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
    M(i,:,:) = K_dB;
    
    % % Make a 2D visualization of all of the recorded signals
    % figure;
    % imagesc(x, t, abs(squeeze(K_dB)).');
    % colorbar;
    % xlabel('lateral distance [mm]');
    % ylabel('time [us]');
    % title('emission #',i)
    % K = K + K_i; % sum up data
    % 
    % K_dB_sum = 20*log10(abs(K).^2/max(abs(K(:))).^2);
    % 
    % % Make a 2D visualization of all of the recorded signals
    % figure;
    % imagesc(x, t, abs(squeeze(K_dB_sum)).');
    % colorbar;
    % xlabel('lateral distance [mm]');
    % ylabel('time [us]');
    % title('compound',i)
end
%%

% Make a 2D visualization of all of the recorded signals
figure;
imagesc(x, t, squeeze(median(M(1:84,:,:))));
colorbar;
xlabel('lateral distance [mm]');
ylabel('time [us]');
title('emission #',i)