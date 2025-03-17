% had a first go and was successful in re-producing figures from the paper;
% Results has spikes with spikeIDs in first column and time-stamps in
% second column; 
% time is: time=0:1/Results.params.RecordingSettings.sampleRate:(Results.params.RecordingSettings.maxRecTime/1000);
% do this to find neurons that respond to the stimulation: spikeIDs =  Results.spikes(Results.spikes(:,2)>Results.params.TissueParams.StimulationOn & ...
%    Results.spikes(:,2)<Results.params.TissueParams.StimulationOff,1)
%% load data from paper
addpath('I:\scripts\Matlab\Vertex_Simulations');
addpath(genpath('I:\scripts\Matlab\tbx\Vertex_git-master'));
F1000_data_dir = 'I:\WT_resub\Vertex_simulations\F1000_Data';
%singlePulseDir = [F1000_data_dir filesep 'pairedpulse\pairedpulse_1001'];
singlePulseDir = [F1000_data_dir filesep '\pairedpulse\cont1sec_1microAmps_1001'];

%singlePulseDir = [F1000_data_dir filesep 'singlepulse_stimon1500_1001'];
%singlePulseDir = [F1000_data_dir filesep 'singlepulse\singlepulse_1001'];

RecordingSettings.saveDir = singlePulseDir;
Results = loadResults(RecordingSettings.saveDir,1);

%% Identify recruited neurons based on increased spiking 100 ms after stimulation
% get spike times
spikeTs = Results.spikes;
neuronIDs = unique(spikeTs(:,1));
% Identify neurons which respond to stimulation by comparing FR 100
% ms after stimulation with baseline FR
bnstp=100;
bins=0:bnstp:2000;
threshld = 5; % this is the threshold for defining a reponsive neuron; it means it needs to increase it's firing rate by that ratio; i.e. 2 means a two-fold increase in FR is required to define this neuron as being responsive
cnt1=0;cnt2=0;
rec_ids=[];notrec_ids=[];
h = waitbar(0, 'Processing...');
for n=1:length(neuronIDs)
    waitbar(n / length(neuronIDs), h, sprintf('Progress: %d%%', round(n / length(neuronIDs) * 100)));
    tmpidx=find(spikeTs(:,1) == neuronIDs(n));
    spktr=hist(spikeTs(tmpidx,2),bins);
    fr_bl=mean(spktr(3:5))./(bnstp/1000);
    if fr_bl == 0
        fr_bl = 0.1;
    end
    fr_stim=mean(spktr(6:16))./(bnstp/1000);
    if (fr_stim / fr_bl) > threshld
        cnt1=cnt1+1;
        rec_ids(cnt1,1)=neuronIDs(n);
    else
        cnt2=cnt2+1;
        notrec_ids(cnt2,1)=neuronIDs(n);
    end
end
close(h)

spikeIDs = rec_ids;

%% Plot electrode location and e-field

setUpBipolarElectrodeStimulation

loadSimulationSettings

figure;
pdeplot3D(model, 'ColorMap', TissueParams.StimulationField.NodalSolution, 'FaceAlpha', 0.8)
hold on;


%% PLot number of Neurons recruited by pulse

% spikeIDs =  Results.spikes(Results.spikes(:,2)>Results.params.TissueParams.StimulationOn(1) & ...
%     Results.params.TissueParams.StimulationOff(1),1);
% spikeTs = Results.spikes;
% rec_ids=[];notrec_ids=[];
% get neuron positions and create 2D histogramm for x and z dimension
% representing the number of neurons per quadrant
all_pos=Results.params.TissueParams.somaPositionMat(:,[1 3]);
x_edges = 0:10:2000;
y_edges = 0:10:2082;
spikeIDs = unique(spikeIDs);

[n_counts, x_bin_edges, y_bin_edges] = histcounts2(all_pos(:,1), all_pos(:,2), x_edges, y_edges);
[a_counts, x_bin_edges, y_bin_edges] = histcounts2(all_pos(spikeIDs,1), all_pos(spikeIDs,2), x_edges, y_edges);

norm_counts = a_counts ./ n_counts;
nanidx=find(isnan(norm_counts));
norm_counts(nanidx)=0;

% Coordinates for the circle
x_circle = 789;
y_circle = 1227;
radius = 10; % Adjust the radius as needed

% Plot the 2D histogram
figure; 
subplot(2,2,1);
imagesc(x_bin_edges, y_bin_edges, n_counts');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist all neurons');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

subplot(2,2,2);
imagesc(x_bin_edges, y_bin_edges, a_counts');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist recruited neurons');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

subplot(2,2,3);
imagesc(x_bin_edges, y_bin_edges, norm_counts');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title(['2D Hist % rec neurons N=' num2str(numel(spikeIDs))]);
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');
% Add location of stimulation electrode
hold on; % Make sure to hold the plot so you can add annotations

% Draw the circle
theta = linspace(0, 2*pi, 100); % Points for circle
x_circle_points = x_circle + radius*cos(theta);
y_circle_points = y_circle + radius*sin(theta);
plot(x_circle_points, y_circle_points, 'r', 'LineWidth', 2); % Red circle

% Apply Gaussian smoothing
sigma = 2; % Standard deviation of the Gaussian filter
snorm_counts = imgaussfilt(norm_counts, sigma);
subplot(2,2,4);
imagesc(x_bin_edges, y_bin_edges, snorm_counts');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist % recruited neurons smoothed');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

% Add location of stimulation electrode
hold on; % Make sure to hold the plot so you can add annotations

% Draw the circle
theta = linspace(0, 2*pi, 100); % Points for circle
x_circle_points = x_circle + radius*cos(theta);
y_circle_points = y_circle + radius*sin(theta);
plot(x_circle_points, y_circle_points, 'r', 'LineWidth', 2); % Red circle

% Add a label for clarity
%text(x_circle, y_circle + 20, 'Electrode', 'Color', 'red', 'FontSize', 12);

%% Plot Spike Histogramms for recruited neurons

bnstp=50;
fr_conv=1000/bnstp;%need this number to convert to firing rate (Hz)
bins=0:bnstp:2000;

spk_matr_rec=zeros(numel(spikeIDs),numel(bins));
h = waitbar(0, 'Processing...');
for n=1:length(spikeIDs)
    waitbar(n / length(spikeIDs), h, sprintf('Progress: %d%%', round(n / length(spikeIDs) * 100)));
    tmpidx=find(spikeTs(:,1) == spikeIDs(n));
    spk_matr_rec(n,:)=hist(spikeTs(tmpidx,2),bins).*20;
end
close(h)

figure;plot(bins,mean(spk_matr_rec,1));xlim([200 2000]);xticks([200 500 1000 1500 2000]);
%% Plot Spike Density of recruited and non-recruited neurons

bnstp=1;
bins=0:1:2000;

spk_matr_rec=zeros(numel(spikeIDs),numel(bins));
spk_dens_rec=spk_matr_rec;

% Define the parameters for the Gaussian window
sigma = 2; % Standard deviation of the Gaussian window
ws = ceil(6*sigma); % Determine the window size based on sigma
gw = normpdf(-ws:ws, 0, sigma); % Generate the Gaussian window
h = waitbar(0, 'Processing...');
for n=1:length(spikeIDs)
    waitbar(n / length(spikeIDs), h, sprintf('Progress: %d%%', round(n / length(spikeIDs) * 100)));
    tmpidx=find(spikeTs(:,1) == spikeIDs(n));
    spk_matr_rec(n,:)=hist(spikeTs(tmpidx,2),bins);
    spk_dens_rec(n,:)= conv(spk_matr_rec(n,:)', gw, 'same')';
end
close(h)

% normalise spike rate to plot % change from baseline
bl=[200 499];
blidx=find(bins == bl(1));
blidx(2)=find(bins == bl(2));
mnfr=mean(spk_dens_rec(:,blidx(1):blidx(2)),2)+0.01;
mnfrmat=repmat(mnfr,[1 size(spk_dens_rec,2)]);
spk_dens_rec_norm=((spk_dens_rec - mnfrmat)./mnfrmat).*100;

% now the same for non-recruited neurons
spk_matr_nrec=zeros(numel(notrec_ids),numel(bins));
spk_dens_nrec=spk_matr_nrec;
h = waitbar(0, 'Processing...');
for n=1:length(notrec_ids)
    waitbar(n / length(notrec_ids), h, sprintf('Progress: %d%%', round(n / length(notrec_ids) * 100)));
    tmpidx=find(spikeTs(:,1) == notrec_ids(n));
    spk_matr_nrec(n,:)=hist(spikeTs(tmpidx,2),bins);
    spk_dens_nrec(n,:)= conv(spk_matr_nrec(n,:)', gw, 'same')';
end
close(h)
% ... normalize
mnfr=mean(spk_dens_nrec(:,blidx(1):blidx(2)),2)+0.01;
mnfrmat=repmat(mnfr,[1 size(spk_dens_nrec,2)]);
spk_dens_nrec_norm=((spk_dens_nrec - mnfrmat)./mnfrmat).*100;

% calculate stanard error around mean
se_rec=nanstd(spk_dens_rec_norm,0,1)./sqrt(numel(spikeIDs));
mn_rec=nanmean(spk_dens_rec_norm,1)+100;
upper=mn_rec+se_rec;
lower=mn_rec-se_rec;
color=[1 0 0];
edge=[1 0 0];
figure;
jbfill(bins,upper,lower,color,edge,0,0.5);
hold on
plot(bins,mn_rec);

se_nrec=std(spk_dens_nrec_norm,0,1)./sqrt(numel(notrec_ids));
mn_nrec=mean(spk_dens_nrec_norm,1)+100;
upper=mn_nrec+se_nrec;
lower=mn_nrec-se_nrec;
color=[0 0 0];
edge=[0 0 0];
jbfill(bins,upper,lower,color,edge,0,0.5);
plot(bins,mn_nrec);
title('Spike Density of recruited and not recruited neurons');
hold off
xlim([100 2000]);
ylim([-20 500]);
% select recruited neurons
%spikeTs_rec = spikeTs()

%% Rasterplot
figure;
% Extract neuron IDs and spike times
neuron_ids = spikeTs(:, 1);
spike_times = spikeTs(:, 2);

% Plot the raster plot
hold on
scatter(spike_times, neuron_ids, 10, 'k', 'filled');
xlabel('Time (ms)');
ylabel('Neuron ID');
title('Raster Plot');
ylim([min(neuron_ids)-1, max(neuron_ids)+1]); % Adjust y-axis limits

% Invert y-axis to display neuron IDs from bottom to top
set(gca, 'YDir', 'reverse');
xlim([0 2000]);

% Identify spikes corresponding to neurons in 'spikeID'
red_indices = ismember(neuron_ids, spikeIDs);

% Plot spikes of neurons in 'spikeID' in red
scatter(spike_times(red_indices), neuron_ids(red_indices), 10, 'k', 'filled'); % Plot spikes of 'spikeID' neurons in red

xlabel('Time (ms)');
ylabel('Neuron ID');
title('Raster Plot');
ylim([min(neuron_ids)-1, max(neuron_ids)+1]); % Adjust y-axis limits
set(gca, 'YDir', 'reverse'); % Invert y-axis


