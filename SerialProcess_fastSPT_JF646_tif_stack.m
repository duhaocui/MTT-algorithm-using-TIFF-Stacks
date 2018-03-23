%   Commit Test
%   SerialProcess_fastSPT_JF646.m
%   Anders Sejr Hansen, August 2016
%%% Edit of Ander Sejr Hansen's code to no longer take nd2 files, but
%%% instead read in Tiff stacks into a 3D matrix and input that into the
%%% MTT algorithm

clear; clc; close all; clearvars -global

%   DESCRIPTION
%   This script takes as input a folder with Tiff stack files and then outputs
%   workspaces with tracked single molecules. It reads the Tiff stack into a 3D
%   matrix. Next, the script feeds 3D matrix into the localization part of the 
%   MTT algorithm (Part 1) and subsequently, the tracked particles are fed 
%   into the tracking part of the MTT algorithm (Part 2). 

%%%%%%%%%%%%%%%%%%%% DEFINE INPUT AND OUTPUT PATHS %%%%%%%%%%%%%%%%%%%%%%%%
% Acuire working directory path
workdir = cd;
% Specify input path with Tiff stack files: The path written as one of the
input_path = uigetdir('C:\Users\Bewersdorf\Desktop\Lukas Fuentes\Data', 'Select Input Data Folder');
% Add path to input Tiff Stacks
addpath(input_path);
% Specify output path for results
output_path= uigetdir('C:\Users\Bewersdorf\Desktop\Lukas Fuentes\Data', 'Select Output Folder');
% uigetdir functions parameters needs to be changed to be specific to your
% system for both input and output paths. This path is just the starting
% path for the uigetdir function.

%%%%% make output folder if it does not already exist
if exist(output_path) == 7
    % OK, output_path exists and is a directory (== 7). 
    disp('The given output folder exists. MATLAB workspaces will be save to:');
    disp(output_path);
else
    mkdir(output_path);
    disp('The given output folder did not exist, but was just created. MATLAB workspaces will be save to:');
    disp(output_path);
end

%%% Parameters of particle tracking %%%
LocalizationError = -6.25; % Localization Error: -6 = 10^-6
EmissionWavelength = 664; % wavelength in nm; consider emission max and filter cutoff
ExposureTime = 4; % in milliseconds
NumDeflationLoops = 0; % Generaly keep this to 0; if you need deflation loops, you are imaging at too high a density;
MaxExpectedD = 10; % The maximal expected diffusion constant for tracking in units of um^2/s;
NumGapsAllowed = 1; % the number of gaps allowed in trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DEFINE STRUCTURED ARRAY WITH ALL THE SPECIFIC SETTINGS FOR LOC AND TRACK
% imaging parameters
impars.PixelSize=0.103; % um per pixel
impars.psf_scale=1.3; % PSF scaling
impars.wvlnth= EmissionWavelength/1000; %emission wavelength in um
impars.NA=1.49; % NA of detection objective
impars.psfStd= impars.psf_scale*0.55*(impars.wvlnth)/impars.NA/1.17/impars.PixelSize/2; % PSF standard deviation in pixels
impars.FrameRate= ExposureTime/1000; %secs
impars.FrameSize= ExposureTime/1000; %secs

% localization parameters
locpars.wn=9; %detection box in pixels
locpars.errorRate= LocalizationError; % error rate (10^-)
locpars.dfltnLoops= NumDeflationLoops; % number of deflation loops
locpars.minInt=0; %minimum intensity in counts
locpars.maxOptimIter= 50; % max number of iterations
locpars.termTol= -2; % termination tolerance
locpars.isRadiusTol=false; % use radius tolerance
locpars.radiusTol=50; % radius tolerance in percent
locpars.posTol= 1.5;%max position refinement
locpars.optim = [locpars.maxOptimIter,locpars.termTol,locpars.isRadiusTol,locpars.radiusTol,locpars.posTol];
locpars.isThreshLocPrec = false;
locpars.minLoc = 0;
locpars.maxLoc = inf;
locpars.isThreshSNR = false;
locpars.minSNR = 0;
locpars.maxSNR = inf;
locpars.isThreshDensity = false;

% tracking parameters
trackpars.trackStart=1;
trackpars.trackEnd=inf;
trackpars.Dmax= MaxExpectedD;
trackpars.searchExpFac=1.2;
trackpars.statWin=10;
trackpars.maxComp=5;
trackpars.maxOffTime=NumGapsAllowed;
trackpars.intLawWeight=0.9;
trackpars.diffLawWeight=0.5;


% add the required functions to the path: needs to be changed to be
% specific to your system.
addpath('C:\Users\Bewersdorf\Desktop\Lukas Fuentes\Single particle tracking\20170727_Maxime\SLIMFAST_batch_fordist');
addpath('C:\Users\Bewersdorf\Desktop\Lukas Fuentes\Single particle tracking\20170727_Maxime\SLIMFAST_batch_fordist\bfmatlab');
disp('added paths for MTT algorithm mechanics, bioformats...');

files = dir(input_path);
Filenames = ''; %for saving the actual file name

for iter = 3:length(files)
    Filenames{iter} = files(iter).name(1:end-4);
end

for iter = 3:length(Filenames)
    disp('-----------------------------------------------------------------');
    tic; 
    disp(['reading in Tiff stack file ', num2str(iter - 2), ' of ', num2str(length(Filenames)-2), ' total Tiff stack files']);

    %%% For uploading Tiff stack %%%
    FileTif = files(iter).name; %uigetfile
    InfoImage = imfinfo(FileTif);
    mImage = InfoImage(1).Width;
    nImage = InfoImage(1).Height;
    NumberImages = length(InfoImage);
    tif_files = zeros(nImage,mImage,NumberImages,'uint16');
    TifLink = Tiff(FileTif, 'r');
    for i = 1:NumberImages
        TifLink.setDirectory(i);
        tif_files(:,:,i) = TifLink.read();
    end
    TifLink.close();
    
    imgs_3d_matrix = tif_files;
    %convert to double, re-scale to 16 bit:
    imgs_3d_double = (2^16-1)*im2double(imgs_3d_matrix);
    
    toc;
    clear img_stack_cell_array cell_array_2d imgs_only_cell_array 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%% MTT ALGORITHM PART 1: LOCALIZE ALL PARTICLES %%%%%%%%%%%%%%
    disp('MTT ALGORITHM PART 1: localize particles in all of the workspaces');
    tic; 
    disp(['localizing all particles in movie number ', num2str(iter-2), ' of ', num2str(length(Filenames)-2)]);
    impars.name = Filenames{iter};
    data = localizeParticles_ASH(input_path,impars, locpars, imgs_3d_matrix);
    %data = localizeParticles_ASH(input_path,impars, locpars, imgs_3d_double);
    toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%% MTT ALGORITHM PART 2: TRACK PARTICLES BETWEEN FRAMES %%%%%%%%%
    disp('MTT ALGORITHM PART 2: track particles between frames');
    tic; 
    disp(['tracking all localized particles from movie ', num2str(iter-2), ' of ', num2str(length(Filenames)-2)]);
    data=buildTracks2_ASH(input_path, data,impars, locpars, trackpars, data.ctrsN, imgs_3d_double);
    toc;
    
    %%%%%%%% SAVE THE TRAJECTORIES TO YOUR STRUCTURED ARRAY FORMAT %%%%%%%%
    tic;
    disp(['saving MATLAB workspace for movie ', num2str(iter-2), ' of ', num2str(length(Filenames)-2)]);
    data_cell_array = data.tr;
    % save meta-data
    settings.Delay = impars.FrameRate;
    settings.px2micron = impars.PixelSize;
    settings.TrackingOptions = trackpars;
    settings.LocOptions = locpars;
    settings.AcquisitionOptions = impars;
    settings.Filename = impars.name;
    settings.Width = size(imgs_3d_matrix,2);
    settings.Height = size(imgs_3d_matrix,1);
    settings.Frames = size(imgs_3d_matrix,3);
    trackedPar = struct;
    for i=1:length(data_cell_array)
        %convert to um:
        trackedPar(1,i).xy =  impars.PixelSize .* data_cell_array{i}(:,1:2);
        trackedPar(i).Frame = data_cell_array{i}(:,3);
        trackedPar(i).TimeStamp = impars.FrameRate.* data_cell_array{i}(:,3);
    end
    disp(['Localized and tracked ', num2str(length(trackedPar)), ' trajectories']);
    % Next four lines are for fixing output file problem (outputs to wrong file) 3/2/18 - LAF
    % mkdir output_path Results_{iter}
    % eval(['mkdir ' output_path ' Results_' num2str(iter);])
    % save_name = fullfile(output_path, Filenames{iter}, '_Tracked.mat');
    % save(save_name, 'trackedPar', 'settings');
    cd(output_path);
    save([Filenames{iter}, '_Tracked.mat'], 'trackedPar', 'settings');
    cd(workdir);
    toc;
    clear imgs_3d_matrix imgs_3d_double data_cell_array trackedPar save_name
    disp('-----------------------------------------------------------------');
end