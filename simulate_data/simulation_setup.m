% run simulations here
% Megumi Hirose (17 June 2025)
% Stores data in a file
% Stored data: 
% 'K_US': Ultrasound data
% 'K_PA': Photoacoustic data
% 'array': transducer array
% 'medium_US': medium used for US simu. Includes reflective target
% 'target' : contains parameters of the reflective target
% 'pulse' : contains parameters of pulse
% 'pml' : perfect matching layer
% 'medium' : medium without the reflective target (use for PA simulation)
% 'kgrid' 


clear variables
clear global
close all hidden
clc

folder.base = pwd;
addpath([folder.base '/subfunctions'])
addpath([folder.base '/output'])

set(0,'DefaultFigureWindowStyle','normal')  % just a preference - you can change this if you like

% open GUI
kgrid_gui_ver12
