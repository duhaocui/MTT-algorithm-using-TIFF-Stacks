%%% This script reads MatLab structures organized in "Ander's Format". That
%%% is, it is a structure where each row corresponds to one trajectory and
%%% column 1 are x-y positions, column 2 are Frame number, and column 3 are
%%% time stamps. Each column contains a Nx2, Nx1, and Nx1, respectively,
%%% matrix, where N is the number of detections for that particle.

%{
Inputs:
- Particle Trajectories in Ander's Format (described above)

Outputs:
- A figure for each trackedPar structure input (presumably equal to the
number of cells you imaged) with all the trajectories of particles tracked
for that trackedPar structure.
%}

% Specify input path for trajectory structures
input_path2 = uigetdir('C:\Users\Bewersdorf\Desktop\Lukas Fuentes\Data', 'Select Input Data Folder');
% Acquire working directory
workdir2 = cd;
% Add input file path
addpath(input_path2);
% Load files from input path
files2 = dir(input_path2);

% Load trajectory structures into cell array
% 'a' starts at 3 to skip the first two columns in 'files' which are
% essentially empty and unavoidable (some weird Windows thing?)
Data = {};
counter = 1;
for a = 3:length(files2)
    if files2(a).name(end-2:end) == 'mat'; % Only loads .mat files
        Data{counter,1} = load(files2(a).name);
        counter = counter + 1;
    end
end

% Plot trajectories
Dleng = (1:length(Data)); % Vector used for title and figure names
plot_title = 'Trajectories in cell #%d';
% Loop over the number of cells/trackPar structures (b), then loop over
% number of trajectories in each trackPar structure (c) and plot each
% trajectory. One figure for each cell will output.
for b = 1:length(Data);
    figure(b);hold on;
    for c = 1:length(Data{b,1}.trackedPar);
        if length(Data{b,1}.trackedPar(c).xy) >= 50; % Adjust threshold to only plot trajectories longer than specified length (in frames)
            plot(Data{b,1}.trackedPar(c).xy(:,1), Data{b,1}.trackedPar(c).xy(:,2));
        end
    end
    % Label axes and title, then save plot in same directory as input files
    xlabel('x-position (\mum)');
    ylabel('y-position (\mum)');
    title(sprintf(plot_title, Dleng(b)));
    cd(input_path2);
    saveas(gcf, sprintf(plot_title, Dleng(b)));
    cd(workdir2);
    hold off;
end