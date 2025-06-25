% run stuff here
% 17 June 2025

clear variables
clear global
% close all hidden
clc

folder.base = pwd;
addpath([folder.base '/subfunctions_ab_corrections']);
addpath([folder.base '/data/Ruu']);
addpath([folder.base '/data/Rxx']);
addpath([folder.base '/GUI/output']);
% 

% %% data with no aberration
% data_file = 'K_US_PA_no_ab.mat'; % contains medium properties and simulated data Ruu
% rxx_file = 'Rxx_no_ab.mat'; % file to store Rxx matrix created from simulated data

% %% data with aberration
data_file = 'K_US_PA2.mat'; % contains medium properties and simulated data Ruu
rxx_file = 'Rxx2.mat'; % file to store Rxx matrix created from simulated data

% if not done yet, run filter_and_Rxx(data_file, rxx_file) to create Rxx
% filter_and_Rxx(data_file, rxx_file);

%%
load(data_file);  % contains Ruu and medium_properties
load(rxx_file);   % contains Rxx

% Image parameters
max_z = kgrid.Nx * kgrid.dx; % [m]
max_x = kgrid.Ny * kgrid.dy; % [m]
zpixels = (-max_z/2:0.05e-3:max_z/2); % [m]
xpixels = (-max_x/2:0.05e-3:max_x/2); % [m]

rx = target.dx/1000;
x_target = (rx + kgrid.x_vec(1) + abs(kgrid.x_vec(1) - array.element.depth(1)) + array.element.height/2)*1000; %[mm]

% image properties
Img.xvec = xpixels;  % [m]
Img.zvec = zpixels;  % [m]
Img.Nx = length(Img.xvec);
Img.Nz = length(Img.zvec);
Img.dx = abs(Img.xvec(1) - Img.xvec(2));
Img.dz = abs(Img.zvec(1) - Img.zvec(2));
Img.c = medium_US.sound_speed_ref;
Img.c_ref = medium_US.sound_speed_ref;
Img.plotflag = 0;
Img.zplot_ind = round((x_target - 1.6)/1000/kgrid.dx); % 1.6 adjustment to account for the aberration - the target is shifted upwards


%%%%%%%%%%% setting up geometry and propagation matrices (G0, T0) %%%%%%%%%

lambda = pulse.wavelength_ref;  % wavelength of emitted pulse at central frequency [m]
% this is equivalent to lambda = medium.sound_speed_ref/pulse.tone_burst_freq_HF;

k = 2*pi/lambda; % wave number [1/m]
theta = -20:0.5:20; % angles of plane-waves [deg]
kx = k*sind(theta).'; % [1/m] transverse wave-number at the central frequency of the emitted pulse
% Willam thesis pg. 71

% T0(kx, x)
T0 = exp(1i*kx*Img.xvec).'/length(kx); % Fourier transform operator : links plane-wave with wave number kx to coordinate on focal plane x
% William thesis, Eq. II.11, pg. 71

%%%%%%%%%%% propagation between focal planes (e.g. u and x)  %%%%%%%%%%%
uvec = array.element.lateral; % [m]
Nu = length(uvec);

[U,X] = meshgrid(uvec.',Img.xvec);

%%%%%% project Rxx onto dual basis %%%%%
Rkk = zeros(length(kx), length(kx),Img.Nz);
Rkx = zeros(Img.Nx,length(kx),Img.Nz); % Rkx(x_in,k_out,z)
Dkx = zeros(size(Rkx));

% we have:  Rxx(x_in,x_out,z)
for ii = 1:Img.Nz
    Rkk(:,:,ii) = T0'*squeeze(Rxx(:,:,ii))*T0;
    Rkx(:,:,ii) = squeeze(Rxx(:,:,ii))*T0;
    Dkx(:,:,ii) = squeeze(Rkx(:,:,ii)).*conj(T0);
end

%% begin reverbaration removal
dk = 1/abs(X(1) - X(end));
[Kin, Kout] = meshgrid(kx, kx);
% Kout = columns kx'
% Kin = rows kx
Dk = abs(Kin + Kout);
Dk_gt = Dk> dk;  % greater than
Dk_lt = Dk <= dk; % less than

% mean abs(Rkk) in depth
R_abs_mean = mean(abs(Rkk), 3);  % size [kk,kk]

gk_sum = sum(R_abs_mean(Dk_gt));
lk_sum = sum(R_abs_mean(Dk_lt));

num_gk = nnz(Dk_gt); % number of greater thans
num_lk = nnz(Dk_lt); % number of less thans

% compute averages
average_gk = gk_sum / num_gk;
average_lk = lk_sum / num_lk;

% compute alpha
alpha = average_gk / average_lk - 1;

% apply reverberation filter
% we have:  Rkk(k_in,k_out,z)
for kin = 1:length(kx)
    for kout = 1:length(kx)
        Rkk(kin, kout, :) = Rkk(kin, kout,:).*(1-alpha*exp(-abs(kx(kout) + kx(kin))^2/dk^2));
    end
end

Ckk = zeros(length(kx), length(kx), Img.Nz);
N = Img.Nx*Img.Nz;  % total number of input focussing points r = (x_in,z)

Rxx_reverb_filtered = Rxx(:,:,:);
for ii = 1:Img.Nz
    Rxx_reverb_filtered(:,:,ii) =T0*Rkk(:,:,ii)*T0';
    Rkx(:,:,ii) = Rxx_reverb_filtered(:,:,ii)*T0;
    Rkx(:,:,ii) =Rkx(:,:,ii)./abs(Rkx(:,:,ii));
    Dkx(:,:,ii) = conj(T0).*Rkx(:,:,ii); % Lambert III.11 pg 105
    % Ckk(:,:,ii) = squeeze(Dkx(:,:,ii))'*squeeze(Dkx(:,:,ii))/N; % C(k_in,k_out)  Eq. III.19 pg 109 William thesis
    % Ckk(:,:,ii) = conj(squeeze(Dkx(:,:,ii)))'*squeeze(Dkx(:,:,ii))/N;

    Ckk(:,:,ii) = conj(Dkx(:,:,ii))'*Dkx(:,:,ii)/N;
    Ckk(:,:,ii) =  Ckk(:,:,ii)./abs(Ckk(:,:,ii));
end

%% SVD filter
Tp = zeros(Img.Nz,Img.Nx, length(kx));
% repeat with every depth z
for ii = 1:Img.Nz
    [U, ~, ~] = svd(Ckk(:,:,ii));
    U1 = U(:,1);             % dominant singular vector
    U1 = U1 ./ abs(U1);      % normalise

    % Tp(ii,:,:) = T0 * diag(U1);  % T0: [Nx × Nk], Tp: [Nz × Nx × Nk]

    Tp(ii,:,:) = T0 .* U1';  % Lambert III.26

end

Rxx_svd_filtered = Rxx_reverb_filtered(:,:,:);
for ii = 1:Img.Nz
    Tpi = squeeze(Tp(ii,:,:));  % [Nx × Nk]
    Tpi = Tpi ./ abs(Tpi); 
    % Tpi = Tpi ./ vecnorm(Tpi);  % optional column-wise normalization
    Rxx_svd_filtered(:,:,ii) = Tpi * Rkk(:,:,ii) * Tpi';
        % Rxx_svd_filtered(:,:,ii) = conj(Tpi) * Rkk(:,:,ii) * conj(Tpi)';
end
% Tp = Up o T0 (Lambert III.26)
% Rxx = transpose(Tp*) x Rkk x Tp* (Lambert III.27)

%%
plot_Rxx(data_file, Rxx, Img, 'Rxx original');
plot_Rxx(data_file, Rxx_reverb_filtered,Img,  'Rxx reverb filtered')
plot_Rxx(data_file, Rxx_svd_filtered, Img, 'Rxx SVD filtered')

%% plotting
plot_slice(Rxx_svd_filtered, T0, Img)





