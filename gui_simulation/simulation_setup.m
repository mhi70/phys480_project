%8

clear variables
clear global
close all hidden
clc

set(0,'DefaultFigureWindowStyle','normal')  % just a preference - you can change this if you like

% set path
directory_content = dir; % contains everything of the current directory
exe_path = directory_content(1).folder; % returns the path that is currently open
addpath(string(exe_path)+''+'/subfunctions')
addpath(string(exe_path)+''+'/output')

% set up simulation environment
kgrid_gui_ver10
