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
clc;clear;
% Specify input path for trajectory structures
input_path2 = uigetdir('C:\Users\Bewersdorf\Desktop\Lukas Fuentes\Data', 'Select Input Data Folder');
% Acquire working directory
workdir2 = cd;
% Add input file path
addpath(input_path2);
% Load files from input path
files2 = dir(input_path2);

%%% Decide Parameters %%%
L = 105; % Minimum number of frames (trajectories below this threshold will be ignored; must be greater than 1)
alpha_bin = 20; % Bin size for alpha histogram
diff_bin = 50; % Bin size for alpha histogram

%%% Load trajectory structures into cell array %%%
% 'a' starts at 3 to skip the first two columns in 'files' which are
% essentially empty and unavoidable (some weird Windows thing?)
Data = {};
counter = 1;
for a = 3:length(files2)
    if files2(a).name(end-2:end) == 'mat' % Only loads .mat files
        Data{counter,1} = load(files2(a).name);
        counter = counter + 1;
    end
end

%% Plot trajectories
Dleng = (1:length(Data)); % Vector used for title and figure names
plot_title = 'Trajectories in cell #%d';
% Loop over the number of cells/trackPar structures (b), then loop over
% number of trajectories in each trackPar structure (c) and plot each
% trajectory. One figure for each cell will output.
for b = 1:length(Data)
    figure(b);hold on;
    for c = 1:length(Data{b,1}.trackedPar)
        if length(Data{b,1}.trackedPar(c).xy) >= L % Adjust threshold to only plot trajectories longer than specified length (in frames)
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

%% Plot MSD and calculate Anomolous Diffusion Coefficient and alpha (power) from MSD ~ Dt^(alpha)

DL = length(Data);
for d = 1:DL
    data_msd = {};
    data_msd = {Data{d,1}.trackedPar.TimeStamp; Data{d,1}.trackedPar.xy}';
    MSD = get_msd_v2(data_msd,d,L,alpha_bin,diff_bin,DL);%,input_path2,workdir);
end

%% 
% This section is the same as the get_msd_v2 function with some tweaked variables for testing (ignore if you are just analyzing data;
% this part is only used for testing if you are editing the code).
figure(); hold on;
x = [];
y = [];
time = [];
counter = 1;
[max_size, ~] = max(cellfun('size', data_msd, 1));
MSD_ALL = NaN(max_size(1),length(data_msd));
DiffCo = [];
for a = 1:length(data_msd)
    trace = [];
    x = [];
    y = [];
    time = [];
    %inds_h = traces(:,4) == pd; % Boolean to find indices for specific cell track
    trace = [data_msd{a,1}, data_msd{a,2}(:,1), data_msd{a,2}(:,2)];
    x = trace(:,2);
    y = trace(:,3);
    time = (trace(:,1)-(trace(1,1)-0.02)); % Normalize times to all start at 0.02 seconds
    %dx = [];
    %dy = [];
    if length(time) > 10
        %poop = time;
        msd_plot = [];
        for t = 1:length(time)-1
            ind = t;
            dx = x(1:end-ind) - x(1+ind:end); % first loop; dx & dy are length(time)-1 long, and decrease by 1 each loop until both are 1 long
            dy = y(1:end-ind) - y(1+ind:end);
            msdi = (mean((dx.^2)+(dy.^2))); %Used to be the mean of this but there is no difference
            MSD_ALL(t,counter) = msdi;
        end
        meanmsd = nanmean(MSD_ALL,2);
        nonans = (MSD_ALL(:,counter)) > 0;
        %suminds = sum(inds_h);
        msd_plot = MSD_ALL(nonans,counter);
        plot(log10(time(1:end-1)),log10(msd_plot));
        poly = polyfit((log10(time(1:end-1))),(log10(msd_plot)),1); % originally only considered 1:6 for each plot (linear range)
        alpha(1,counter) = poly(1);
        diff = 10^(poly(2));
        DiffCo(1,counter) = diff;
        counter = counter + 1;
        %disp(a);
    end
end

xlabel('log(time(ms))');
ylabel('log(Squared Displacement(\mum^2))');
title(sprintf('Cell #%d squared displacements', 1));
hold off;
%{
cd(input_path);
saveas(gcf,[input_path '/subFolderName/myFig.fig']);
cd(workdir);
%}
figure(); hold on;
plot_time = linspace(0,0.02,max_size(1))';
plot(log10(plot_time(1:end-1)),log10(meanmsd(1:end-1)));
xlabel('log(time(s))');
ylabel('log(MSD(\mu^2))');
title(sprintf('Cell #%d MSD', 1)); hold off;
%{
figure(); hold on;
histogram(alpha,20);
xlabel('\alpha');
ylabel('Events');
title(sprintf('Cell #%d alpha distribution', d)); hold off;
figure(); hold on;
histogram(DiffCo,20);
xlabel('Diffusion coefficient');
ylabel('Events');
title(sprintf('Cell #%d diffusion coefficient distribution', d)); hold off;
%}
MSD = nanmean(MSD_ALL,2);