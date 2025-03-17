% had a first go and was successful in re-producing figures from the paper;
% Results has spikes with spikeIDs in first column and time-stamps in
% second column; 
% time is: time=0:1/Results.params.RecordingSettings.sampleRate:(Results.params.RecordingSettings.maxRecTime/1000);
% do this to find neurons that respond to the stimulation: spikeIDs =  Results.spikes(Results.spikes(:,2)>Results.params.TissueParams.StimulationOn & ...
%    Results.spikes(:,2)<Results.params.TissueParams.StimulationOff,1)
%% load data from paper
addpath(genpath('I:\scripts\Matlab\tbx\Vertex_git-master'));
F1000_data_dir = 'I:\WT_resub\Vertex_simulations\F1000_Data';
singlePulseDir = [F1000_data_dir filesep '\pairedpulse\cont1sec_1001'];
RecordingSettings.saveDir = singlePulseDir;
Results.exp = loadResults(RecordingSettings.saveDir,1);

singlePulseDir = [F1000_data_dir filesep '\pairedpulse\shrt1sec_1001'];
RecordingSettings.saveDir = singlePulseDir;
Results.bl = loadResults(RecordingSettings.saveDir,1);

%% Identify recruited neurons based on increased spiking comparing exp and baseline (bl)
% get spike times
spikeTs.bl = Results.bl.spikes;
spikeTs.exp = Results.exp.spikes;
neuronIDs_exp = unique(spikeTs.exp(:,1));% if same random seed is used, then neuronIDs will be identical for exp and bl so we only need to grab it once
neuronIDs_bl = unique(spikeTs.bl(:,1));% if same random seed is used, then neuronIDs will be identical for exp and bl so we only need to grab it once
neuronIDs = union(neuronIDs_exp, neuronIDs_bl);

% Identify neurons which respond to stimulation by comparing FR 100
% ms after stimulation with baseline FR
threshld = 10; % this is the threshold for defining a reponsive neuron; it means it needs to increase it's firing rate by that ratio; i.e. 2 means a two-fold increase in FR is required to define this neuron as being responsive
cnt1=0;cnt2=0;cnt3=0;
rec_ids=[];notrec_ids=[];
h = waitbar(0, 'Processing...');
for n=1:length(neuronIDs)
    waitbar(n / length(neuronIDs), h, sprintf('Progress: %d%%', round(n / length(neuronIDs) * 100)));
    fr_bl=numel(find(spikeTs.bl(:,1) == neuronIDs(n)))./2;
    fr_exp=numel(find(spikeTs.exp(:,1) == neuronIDs(n)))./2;
    fr_ratio=fr_exp / fr_bl;
    if fr_ratio > threshld
        cnt1=cnt1+1;
        excited_ids(cnt1,1)=neuronIDs(n);
    elseif fr_ratio < 1
        cnt2=cnt2+1;
        inhib_ids(cnt2,1)=neuronIDs(n);
    else
        cnt3=cnt3+1;
        nonstim_ids(cnt3,1)=neuronIDs(n);
    end
end
close(h)

spikeIDs_e = excited_ids;
spikeIDs_i = inhib_ids;
spikeIDs_n = nonstim_ids;

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

% Apply Gaussian smoothing
sigma = 2; % Standard deviation of the Gaussian filter
all_pos=Results.exp.params.TissueParams.somaPositionMat(:,[1 3]);
x_edges = 0:10:2000;
y_edges = 0:10:2082;
spikeIDs_e = unique(spikeIDs_e);
spikeIDs_i = unique(spikeIDs_i);

[n_counts, x_bin_edges, y_bin_edges] = histcounts2(all_pos(:,1), all_pos(:,2), x_edges, y_edges);
[a_counts_e, x_bin_edges, y_bin_edges] = histcounts2(all_pos(spikeIDs_e,1), all_pos(spikeIDs_e,2), x_edges, y_edges);
[a_counts_i, x_bin_edges, y_bin_edges] = histcounts2(all_pos(spikeIDs_i,1), all_pos(spikeIDs_i,2), x_edges, y_edges);

norm_counts_e = a_counts_e ./ n_counts;
nanidx=find(isnan(norm_counts_e));
norm_counts_e(nanidx)=0;
snorm_counts_e = imgaussfilt(norm_counts_e, sigma);

norm_counts_i = a_counts_i ./ n_counts;
nanidx=find(isnan(norm_counts_i));
norm_counts_i(nanidx)=0;
snorm_counts_i = imgaussfilt(norm_counts_i, sigma);

% Plot the 2D histogram for excited neurons
figure; 
subplot(2,2,1);
imagesc(x_bin_edges, y_bin_edges, norm_counts_e');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist % excited neurons');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

subplot(2,2,3);
imagesc(x_bin_edges, y_bin_edges, norm_counts_i');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist % excited neurons');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

subplot(2,2,2);
imagesc(x_bin_edges, y_bin_edges, snorm_counts_e');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist % inhibited neurons smoothed');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

% Add location of stimulation electrode
hold on; % Make sure to hold the plot so you can add annotations

% Coordinates for the circle
x_circle = 789;
y_circle = 1227;
radius = 10; % Adjust the radius as needed

% Draw the circle
theta = linspace(0, 2*pi, 100); % Points for circle
x_circle_points = x_circle + radius*cos(theta);
y_circle_points = y_circle + radius*sin(theta);
plot(x_circle_points, y_circle_points, 'r', 'LineWidth', 2); % Red circle
hold off

subplot(2,2,4);
imagesc(x_bin_edges, y_bin_edges, snorm_counts_i');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist % inhibited neurons smoothed');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

% Add location of stimulation electrode
hold on; % Make sure to hold the plot so you can add annotations

% Coordinates for the circle
x_circle = 789;
y_circle = 1227;
radius = 10; % Adjust the radius as needed

% Draw the circle
theta = linspace(0, 2*pi, 100); % Points for circle
x_circle_points = x_circle + radius*cos(theta);
y_circle_points = y_circle + radius*sin(theta);
plot(x_circle_points, y_circle_points, 'r', 'LineWidth', 2); % Red circle
hold off

%% Plot FR change bl to exp as 2D histogramm to have better information where spiking increases in response to stimulation
% get spike times
h = waitbar(0, 'Processing...');
for n=1:length(neuronIDs)
    waitbar(n / length(neuronIDs), h, sprintf('Progress: %d%%', round(n / length(neuronIDs) * 100)));
    tmpbl=find(spikeTs.bl(:,1)==neuronIDs(n));
    if ~isempty(tmpbl)
        spikeTs.bl(tmpbl,3:4)=repmat(all_pos(neuronIDs(n),:),[length(tmpbl) 1]);
    end
    tmpexp=find(spikeTs.exp(:,1)==neuronIDs(n));
    if ~isempty(tmpexp)
        spikeTs.exp(tmpexp,3:4)=repmat(all_pos(neuronIDs(n),:),[length(tmpexp) 1]);
    end    
end
close(h)

% Apply Gaussian smoothing
sigma = 2; % Standard deviation of the Gaussian filter
spikes_bl=spikeTs.bl(:,3:4);
spikes_exp=spikeTs.exp(:,3:4);
x_edges = 0:10:2000;
y_edges = 0:10:2082;
[counts_bl, x_bin_edges, y_bin_edges] = histcounts2(spikes_bl(:,1), spikes_bl(:,2), x_edges, y_edges);
[counts_exp, x_bin_edges, y_bin_edges] = histcounts2(spikes_exp(:,1), spikes_exp(:,2), x_edges, y_edges);

% calculate normalised firing rate histos
norm_counts = counts_exp - counts_bl ./ counts_bl;
nanidx=find(isnan(norm_counts));
norm_counts(nanidx)=0;

snorm_counts = imgaussfilt(norm_counts, sigma);
s_counts_bl = imgaussfilt(counts_bl, sigma);
s_counts_exp = imgaussfilt(counts_exp, sigma);

% plot
figure; 
subplot(2,2,1);
imagesc(x_bin_edges, y_bin_edges, counts_bl');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist spikes baseline');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

subplot(2,2,2);
imagesc(x_bin_edges, y_bin_edges, counts_exp');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist % stimulation');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

subplot(2,2,3);
imagesc(x_bin_edges, y_bin_edges, norm_counts', [0 1000]);
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist % stimulation');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

subplot(2,2,4);
imagesc(x_bin_edges, y_bin_edges, snorm_counts', [0 1000]);
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Hist % stimulation');
ax=gca;
set(ax, 'YDir', 'normal');
set(ax, 'XDir', 'reverse');

%% Plot normalized firing rate per 100 ms bin

time_bin =0:100:2000;
figure;
for k=1:numel(time_bin)-1
    sigma = 2; % Standard deviation of the Gaussian filter
    tmp_bl=find(spikeTs.bl(:,2) > time_bin(k) &  spikeTs.bl(:,2) < time_bin(k+1));
    spikes_bl=spikeTs.bl(tmp_bl,3:4);
    tmp_exp=find(spikeTs.exp(:,2) > time_bin(k) &  spikeTs.exp(:,2) < time_bin(k+1));
    spikes_exp=spikeTs.exp(tmp_exp,3:4);
    x_edges = 0:10:2000;
    y_edges = 0:10:2082;
    [counts_bl, x_bin_edges, y_bin_edges] = histcounts2(spikes_bl(:,1), spikes_bl(:,2), x_edges, y_edges);
    [counts_exp, x_bin_edges, y_bin_edges] = histcounts2(spikes_exp(:,1), spikes_exp(:,2), x_edges, y_edges);

    % calculate normalised firing rate histos
    norm_counts = counts_exp - counts_bl ./ counts_bl;
    nanidx=find(isnan(norm_counts));
    norm_counts(nanidx)=0;

    snorm_counts = imgaussfilt(norm_counts, sigma);
    s_counts_bl = imgaussfilt(counts_bl, sigma);
    s_counts_exp = imgaussfilt(counts_exp, sigma);
    
    subplot(4,5,k);
    hold on
    imagesc(x_bin_edges, y_bin_edges, snorm_counts', [0 50]);
    colorbar;
    xlabel('X-axis');
    ylabel('Y-axis');
    title([num2str(time_bin(k)) 'ms']);
    ax=gca;
    set(ax, 'YDir', 'normal');
    set(ax, 'XDir', 'reverse');
    
    % Coordinates for the circle
    x_circle = 789;
    y_circle = 1227;
    radius = 10; % Adjust the radius as needed
    
    % Draw the circle
    theta = linspace(0, 2*pi, 100); % Points for circle
    x_circle_points = x_circle + radius*cos(theta);
    y_circle_points = y_circle + radius*sin(theta);
    plot(x_circle_points, y_circle_points, 'r', 'LineWidth', 2); % Red circle
    ylim([0 2000])
    hold off
    
    
end

%% Plot Spike Density of recruited and non-recruited neurons

bnstp=1;
bins=0:bnstp:2000;

spk_matr_e=zeros(numel(spikeIDs_e),numel(bins));
spk_dens_e=spk_matr_e;

spk_matr_i=zeros(numel(spikeIDs_i),numel(bins));
spk_dens_i=spk_matr_e;

spk_matr_n=zeros(numel(spikeIDs_n),numel(bins));
spk_dens_n=spk_matr_n;

% Define the parameters for the Gaussian window
sigma = 2; % Standard deviation of the Gaussian window
ws = ceil(6*sigma); % Determine the window size based on sigma
gw = normpdf(-ws:ws, 0, sigma); % Generate the Gaussian window

% Calculate spike density for excited neurons
h = waitbar(0, 'Processing...');
for n=1:length(spikeIDs_e)
    waitbar(n / length(spikeIDs_e), h, sprintf('Progress: %d%%', round(n / length(spikeIDs_e) * 100)));
    tmpidx=find(spikeTs.exp(:,1) == spikeIDs_e(n));
    spk_matr_e_exp(n,:)=hist(spikeTs.exp(tmpidx,2),bins);
    tmpidx=find(spikeTs.bl(:,1) == spikeIDs_e(n));
    spk_matr_e_bl(n,:)=hist(spikeTs.bl(tmpidx,2),bins);
    spk_matr_e_diff(n,:)=spk_matr_e_exp(n,:) - spk_matr_e_bl(n,:);
    spk_dens_e_diff(n,:)= conv(spk_matr_e_diff(n,:)', gw, 'same')';
    spk_dens_e_raw(n,:)= conv(spk_matr_e_exp(n,:)', gw, 'same')';
end
close(h)

figure;plot(bins, mean(spk_dens_e_raw,1), 'b');hold on; 
plot(bins, mean(spk_dens_e_diff,1), 'r');hold off; 

% Calculate spike density for inhibited neurons
h = waitbar(0, 'Processing...');
for n=1:length(spikeIDs_i)
    waitbar(n / length(spikeIDs_i), h, sprintf('Progress: %d%%', round(n / length(spikeIDs_i) * 100)));
    tmpidx=find(spikeTs.exp(:,1) == spikeIDs_i(n));
    spk_matr_i_exp(n,:)=hist(spikeTs.exp(tmpidx,2),bins);
    tmpidx=find(spikeTs.bl(:,1) == spikeIDs_i(n));
    spk_matr_i_bl(n,:)=hist(spikeTs.bl(tmpidx,2),bins);
    spk_matr_i_diff(n,:)=spk_matr_i_exp(n,:) - spk_matr_i_bl(n,:);
    spk_dens_i_diff(n,:)= conv(spk_matr_i_diff(n,:)', gw, 'same')';
    spk_dens_i_raw(n,:)= conv(spk_matr_i_exp(n,:)', gw, 'same')';
end
close(h)

figure;plot(bins, mean(spk_dens_i_raw,1), 'b');hold on; 
plot(bins, mean(spk_dens_i_diff,1), 'r');hold off; 

% Calculate spike density for non-stimulated neurons
h = waitbar(0, 'Processing...');
for n=1:length(spikeIDs_n)
    waitbar(n / length(spikeIDs_n), h, sprintf('Progress: %d%%', round(n / length(spikeIDs_n) * 100)));
    tmpidx=find(spikeTs.exp(:,1) == spikeIDs_n(n));
    spk_matr_n_exp(n,:)=hist(spikeTs.exp(tmpidx,2),bins);
    tmpidx=find(spikeTs.bl(:,1) == spikeIDs_n(n));
    spk_matr_n_bl(n,:)=hist(spikeTs.bl(tmpidx,2),bins);
    spk_matr_n_diff(n,:)=spk_matr_n_exp(n,:) - spk_matr_n_bl(n,:);
    spk_dens_n_diff(n,:)= conv(spk_matr_n_diff(n,:)', gw, 'same')';
    spk_dens_n_raw(n,:)= conv(spk_matr_n_exp(n,:)', gw, 'same')';
end
close(h)

figure;plot(bins, mean(spk_dens_n_raw,1), 'b');hold on; 
plot(bins, mean(spk_dens_n_diff,1), 'r');hold off; 

% calculate stanard error around mean
se_e=nanstd(spk_dens_e_diff,0,1)./sqrt(numel(spikeIDs_e));
mn_e=nanmean(spk_dens_e_diff,1);
upper_e=mn_e+se_e;
lower_e=mn_e-se_e;

se_i=nanstd(spk_dens_i_diff,0,1)./sqrt(numel(spikeIDs_i));
mn_i=nanmean(spk_dens_i_diff,1);
upper_i=mn_i+se_i;
lower_i=mn_i-se_i;

se_n=nanstd(spk_dens_n_diff,0,1)./sqrt(numel(spikeIDs_n));
mn_n=nanmean(spk_dens_n_diff,1);
upper_n=mn_n+se_n;
lower_n=mn_n-se_n;

%% Plot
color=[0 0 1];
edge=[0 0 1];
figure;
jbfill(bins,upper_e,lower_e,color,edge,0,0.5);
hold on
plot(bins,mn_e, 'Color', color);

color=[0 0 0];
edge=[0 0 0];
jbfill(bins,upper_n,lower_n,color,edge,0,0.5);
plot(bins,mn_n, 'Color', color);

color=[1 0 0];
edge=[1 0 0];
jbfill(bins,upper_i,lower_i,color,edge,0,0.5);
plot(bins,mn_i, 'Color', color);

title('Spike Density of recruited and not recruited neurons');
hold off
xlim([0 2000]);
xticks([0 500 1000 1500 2000]);
%ylim([-20 500]);
% select recruited neurons
%spikeTs_rec = spikeTs()

%% Rasterplot
figure;
% Extract neuron IDs and spike times
neuron_ids = spikeTs.exp(:, 1);
spike_times = spikeTs.exp(:, 2);

% Plot the raster plot
hold on
scatter(spike_times, neuron_ids, 3, 'k', 'filled');
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
scatter(spike_times(red_indices), neuron_ids(red_indices), 3, 'k', 'filled'); % Plot spikes of 'spikeID' neurons in red

xlabel('Time (ms)');
ylabel('Neuron ID');
title('Raster Plot');
ylim([min(neuron_ids)-1, max(neuron_ids)+1]); % Adjust y-axis limits
set(gca, 'YDir', 'reverse'); % Invert y-axis
xticks([0 500 1000 1500 2000]);

