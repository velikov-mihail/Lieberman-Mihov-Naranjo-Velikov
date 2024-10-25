
clear
clc

% Start with the default path
restoredefaultpath; 

% Path to the MATLAB asset pricing package
matlabPackagePath = 'D:\Published Repos\Lieberman-Mihov-Naranjo-Velikov\'; 

% Path to the code for the paper
paperCodePath = 'D:\Published Repos\Lieberman-Mihov-Naranjo-Velikov\Lieberman-Mihov-Naranjo-Velikov\'; 

% Path to the folder with inputs that should contain 
inputsPath = 'D:\Published Repos\Lieberman-Mihov-Naranjo-Velikov\Lieberman-Mihov-Naranjo-Velikov\Inputs\'; 

% Add the relevant folders (with subfolders) to the path
addpath(genpath([matlabPackagePath, 'Data']))
addpath(genpath([matlabPackagePath, 'Functions']))
addpath(genpath([matlabPackagePath, 'Library Update']))
addpath(genpath([paperCodePath]))
addpath(genpath([inputsPath]))

% Navigate to the paper folder
cd(paperCodePath)

%% Add the directories

% Check if the /Data/ directory exists
if ~exist([pwd,'Data'], 'dir')
    mkdir(['Data'])
end

% Check if the /Results/ directory exists
if ~exist([pwd,'Results'], 'dir')
    mkdir(['Results'])
end

% Check if the /Figures/ directory exists
if ~exist([pwd,'Figures'], 'dir')
    mkdir(['Figures'])
end

% Make sure we add those to the path if we created them
addpath(genpath(pwd));

%% Start a log file

startLogFile(paperCodePath, 'lnmv')

%% Make the data

run('make_data.m');

%% Make the figures

run('make_figures.m');

%% Make the tables

run('make_tables.m');

%% End the log file

diary off
