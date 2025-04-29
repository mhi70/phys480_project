%8

clear variables
clear global
close all hidden
clc

Laura = false;
run_sim = true; 

if Laura
    folder.base = 'C:\Users\lco137\OneDrive - University of Canterbury\Project management\Imaging Team NZ\Megumi honours 2025\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_april15';
    cd(folder.base)
    addpath([folder.base '/subfunctions']);
else
    folder.base = 'C:\Users\mhi70\OneDrive - University of Canterbury\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_april15';
    %folder.base='/Users/megumihirose/Desktop/PHYS480_project/AcousticBeamformingLab_scripts';
    % folder.base = 'C:\Users\mhi70\OneDrive - University of Canterbury\Megumi_photoacoustic_aberration_shared_2025\Beamforming & imaging lab\megumi_code_march_18';
    cd(folder.base)
    addpath([folder.base '/subfunctions']); % or wherever the subfunctions folder is
end
set(0,'DefaultFigureWindowStyle','normal')  % jucst a preference - you can change this if you like

% set up simulation environment
kgrid_gui_ver10
