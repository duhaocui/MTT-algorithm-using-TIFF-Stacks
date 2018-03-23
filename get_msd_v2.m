
function MSD = get_msd_v2(data_msd,n,L,alpha_bin,diff_bin,DL)%input_path,workdir)
figure(DL+n); hold on;
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
    msd_plot = [];
    %inds_h = traces(:,4) == pd; % Boolean to find indices for specific cell track
    if length(data_msd{a,1}) > L
        expo_time = data_msd{a,1}(2)-data_msd{a,1}(1);
        trace = [data_msd{a,1}, data_msd{a,2}(:,1), data_msd{a,2}(:,2)];
        x = trace(:,2);
        y = trace(:,3);
        time = (trace(:,1)-(trace(1,1) - expo_time)); % Normalize times to all start at 0.02 seconds
        dx = [];
        dy = [];
        for t = 1:length(time)-1
            ind = t;
            dx = x(1:end-ind) - x(1+ind:end);
            dy = y(1:end-ind) - y(1+ind:end);
            msdi = (mean((dx.^2)+(dy.^2)));
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
xlabel('log(lag time(s))');
ylabel('log(Squared Displacement(\mum^2))');
title(sprintf('Cell #%d squared displacements', n));
hold off;
%{
cd(input_path);
saveas(gcf,[input_path '/subFolderName/myFig.fig']);
cd(workdir);
%}
figure(2*DL+n); hold on;
plot_time = linspace(expo_time,expo_time*max_size(1),max_size(1))';
plot(log10(plot_time(1:end-1)),log10(meanmsd(1:end-1)));
xlabel('log(lag time(s))');
ylabel('log(MSD(\mum^2))');
title(sprintf('Cell #%d MSD', n));
%{
figure(); hold on;
histogram(alpha,alpha_bin);
xlabel('\alpha');
ylabel('Events');
title(sprintf('Cell #%d alpha distribution', n)); hold off;
figure(); hold on;
histogram(DiffCo,diff_bin);
xlabel('Diffusion coefficient');
ylabel('Events');
title(sprintf('Cell #%d diffusion coefficient distribution', n)); hold off;
%}
MSD = nanmean(MSD_ALL,2);