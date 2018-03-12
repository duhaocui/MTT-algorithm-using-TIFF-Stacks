%%% This script reads MatLab structures organized in "Ander's Format". That
%%% is, it is a structure where each row corresponds to one trajectory and
%%% column 1 are x-y positions, column 2 are Frame number, and column 3 are
%%% time stamps. Each column contains a Nx2, Nx1, and Nx1, respectively,
%%% matrix, where N is the number of detections for that particle.

%{
Inputs:
- Particle Trajectories in Ander's Format (described above)

Outputs:
- Traces of the particle trajectories
%}

% Specify input path for trajectory structures
input_path2 = uigetdir('C:\Users\Bewersdorf\Desktop\Lukas Fuentes\Data', 'Select Input Data Folder');
% Add input file path
addpath(input_path2);
% Load data from input path
files2 = dir(input_path2);
%Filenames2 = ''; %for saving the actual file name

% 'a' starts at 3 to skip the first two columns in 'files' which are
% essentially empty and unavoidable (some weird Windows thing?)
%{
for a = 3:length(files2)
    %Filenames2{a-2} = files2(a).name(1:end);
    s = load(files2(a).name);
    disp(s);
end
%}
s = load(files2(3:4).name);