% function filter_and_Rxx(data_file, output_file)
% Processing data from full matrix capture US simulation. 
%
% Inputs:
%
% data_file : name of file containing simulated data 
% output_file : name of file where output will be stored. 
%
% Function loads data in data_file and applies frequency (HANN) and time (Tukey) filters. 
% Saved variables: 
% save(output_file,'Rxx');
%
% Megumi Hirose (May 2025)
%

function filter_and_Rxx(data_file, output_file)
load(data_file);

%% image params used in US B-mode images

max_z = kgrid.Nx * kgrid.dx; % [m]
max_x = kgrid.Ny * kgrid.dy; % [m]
zpixels = (-max_z/2:0.05e-3:max_z/2); % [m]
xpixels = (-max_x/2:0.05e-3:max_x/2); % [m]

rx = target.dx/1000;
x_target = (rx + kgrid.x_vec(1) + abs(kgrid.x_vec(1) - array.element.depth(1)) + array.element.height/2)*1000; %[mm]

K_US = tukey_hann_filter(data_file);

% populate Img

Img.xvec = xpixels;  % [m]
Img.zvec = zpixels;  % [m]
Img.Nx = length(Img.xvec);
Img.Nz = length(Img.zvec);
Img.c = medium_US.sound_speed_ref;
Img.c_ref = medium_US.sound_speed_ref;
Img.plotflag = 0;
Img.zplot_ind = round((x_target - 1.6)/1000/kgrid.dx); % 1.6 adjustment to account for the aberration - the target is shifted upwards

% create focussed reflection matrix Rxx
[Rxx, ~] = Rxx_DAS_focusing(K_US,kgrid,array,Img);
% gives Rxx(x_in,x_out,z)

save(strcat('data/Rxx/',output_file),'Rxx', 'Img');
end



