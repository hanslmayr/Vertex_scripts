%% Plot to compare results from different stimulation intensities
% set path
addpath('I:\scripts\Matlab\Vertex_Simulations');
addpath(genpath('I:\scripts\Matlab\tbx\Vertex_git-master'));
F1000_data_dir = 'I:\WT_resub\Vertex_simulations\F1000_Data';
intens={'20microAmps';'40microAmps';'60microAmps';'80microAmps'};
for k=1:numel(intens)
    singlePulseDir = [F1000_data_dir filesep '\pairedpulse\cont1sec_' intens{k,1}, '_1001'];
    RecordingSettings.saveDir = singlePulseDir;
    Results(k,1) = loadResults(RecordingSettings.saveDir,1);
end

%% Identify recruited neurons based on increased spiking 100 ms after stimulation
figure; 

for k=1:numel(Results)
    
    % get spike times
    spikeTs{k,1} = Results(k,1).spikes;
    neuronIDs = unique(spikeTs{k,1}(:,1));
    % Identify neurons which respond to stimulation by comparing FR 100
    % ms after stimulation with baseline FR
    bnstp=100;
    bins=0:bnstp:2000;
    threshld = 10; % this is the threshold for defining a reponsive neuron; it means it needs to increase it's firing rate by that ratio; i.e. 2 means a two-fold increase in FR is required to define this neuron as being responsive
    cnt1=0;cnt2=0;
    rec_ids=[];notrec_ids=[];
    h = waitbar(0, 'Processing...');
    for n=1:length(neuronIDs)
        waitbar(n / length(neuronIDs), h, sprintf('Progress: %d%%', round(n / length(neuronIDs) * 100)));
        tmpidx=find(spikeTs{k,1}(:,1) == neuronIDs(n));
        spktr=hist(spikeTs{k,1}(tmpidx,2),bins);
        fr_bl=mean(spktr(1:5))./(bnstp/1000);
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
    
    spikeIDs{k,1} = rec_ids;

%     % Plot electrode location and e-field
% 
%     setUpBipolarElectrodeStimulation
% 
%     loadSimulationSettings
% 
%     figure;
%     pdeplot3D(model, 'ColorMap', TissueParams.StimulationField.NodalSolution, 'FaceAlpha', 0.8)
%     hold on;


    % PLot number of Neurons recruited by pulse

    all_pos=Results(k,1).params.TissueParams.somaPositionMat(:,[1 3]);
    x_edges = 0:10:2000;
    y_edges = 0:10:2082;
    spikeIDs{k,1} = unique(spikeIDs{k,1});

    [n_counts, x_bin_edges, y_bin_edges] = histcounts2(all_pos(:,1), all_pos(:,2), x_edges, y_edges);
    [a_counts, x_bin_edges, y_bin_edges] = histcounts2(all_pos(spikeIDs{k,1},1), all_pos(spikeIDs{k,1},2), x_edges, y_edges);

    norm_counts = a_counts ./ n_counts;
    nanidx=find(isnan(norm_counts));
    norm_counts(nanidx)=0;

    % Coordinates for the circle
    x_circle = 789;
    y_circle = 1227;
    radius = 10; % Adjust the radius as needed

    % Plot the 2D histogram
     % Apply Gaussian smoothing
    sigma = 2; % Standard deviation of the Gaussian filter
    snorm_counts = imgaussfilt(norm_counts, sigma);
    subplot(2,4,k);
    imagesc(x_bin_edges, y_bin_edges, snorm_counts');
    colorbar;
    xlabel('X-axis');
    ylabel('Y-axis');
    title(['% recruited neurons I=' num2str(intens{k,1}(1,1)) 'mA']);
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
    hold off
    N_rec_nrns(k,1)=numel(spikeIDs{k,1});
end

figure;bar(N_rec_nrns);
%% Plot Firing rate per stimulation intensity

bnstp=50;
fr_conv=1000/bnstp;%need this number to convert to firing rate (Hz)
bins=0:bnstp:2000;
figure;
color=[1, 0.6, 0.6;1, 0, 0;0.5, 0, 0; 0.5, 0, 0.5];
for k=1:numel(Results) 
    act_spks=spikeIDs{k,1};
    act_spkTs=spikeTs{k,1};
    spk_matr_rec=zeros(numel(act_spks),numel(bins));
    h = waitbar(0, 'Processing...');
    for n=1:length(act_spks)
        waitbar(n / length(act_spks), h, sprintf('Progress: %d%%', round(n / length(act_spks) * 100)));
        tmpidx=find(act_spkTs(:,1) == act_spks(n));
        spk_matr_rec(n,:)=hist(act_spkTs(tmpidx,2),bins).*20;
    end
    close(h)
    plot(bins,mean(spk_matr_rec,1),'Color', color(k,:));xlim([200 2000]);xticks([200 500 1000 1500 2000]);
    hold on
end



%% Plot Spike Density of recruited and non-recruited neurons
figure;
color=[1, 0.6, 0.6;1, 0, 0;0.5, 0, 0; 0.5, 0, 0.5];
for k=1:numel(Results) 
    
    bnstp=1;
    bins=0:1:2000;
    spk_matr_rec=[];
    spk_dens_rec=[];
    spk_matr_rec=zeros(numel(spikeIDs{k,1}),numel(bins));
    spk_dens_rec=spk_matr_rec;

    % Define the parameters for the Gaussian window
    sigma = 5; % Standard deviation of the Gaussian window
    ws = ceil(6*sigma); % Determine the window size based on sigma
    gw = normpdf(-ws:ws, 0, sigma); % Generate the Gaussian window
    h = waitbar(0, 'Processing...');
    for n=1:length(spikeIDs{k,1})
        waitbar(n / length(spikeIDs{k,1}), h, sprintf('Progress: %d%%', round(n / length(spikeIDs{k,1}) * 100)));
        tmpidx=find(spikeTs(:,1) == spikeIDs{k,1}(n));
        spk_matr_rec(n,:)=hist(spikeTs(tmpidx,2),bins);
        spk_dens_rec(n,:)= conv(spk_matr_rec(n,:)', gw, 'same')';
    end
    close(h)

    % normalise spike rate to plot % change from baseline
    bl=[1 499];
    blidx=find(bins == bl(1));
    blidx(2)=find(bins == bl(2));
    mnfr=mean(spk_dens_rec(:,blidx(1):blidx(2)),2)+0.01;
    mnfrmat=repmat(mnfr,[1 size(spk_dens_rec,2)]);
    spk_dens_rec_norm=((spk_dens_rec - mnfrmat)./mnfrmat).*100;

%     % now the same for non-recruited neurons
%     spk_matr_nrec=zeros(numel(notrec_ids),numel(bins));
%     spk_dens_nrec=spk_matr_nrec;
%     h = waitbar(0, 'Processing...');
%     for n=1:length(notrec_ids)
%         waitbar(n / length(notrec_ids), h, sprintf('Progress: %d%%', round(n / length(notrec_ids) * 100)));
%         tmpidx=find(spikeTs(:,1) == notrec_ids(n));
%         spk_matr_nrec(n,:)=hist(spikeTs(tmpidx,2),bins);
%         spk_dens_nrec(n,:)= conv(spk_matr_nrec(n,:)', gw, 'same')';
%     end
%     close(h)
%     % ... normalize
%     mnfr=mean(spk_dens_nrec(:,blidx(1):blidx(2)),2)+0.01;
%     mnfrmat=repmat(mnfr,[1 size(spk_dens_nrec,2)]);
%     spk_dens_nrec_norm=((spk_dens_nrec - mnfrmat)./mnfrmat).*100;
% 
    % calculate stanard error around mean
    se_rec=nanstd(spk_dens_rec_norm,0,1)./sqrt(numel(spikeIDs{k,1}));
    mn_rec=nanmean(spk_dens_rec_norm,1)+100;
    upper=mn_rec+se_rec;
    lower=mn_rec-se_rec;
    %color=[1 0 0];
    edge=[1 0 0];
    subplot(2,3,3:4);
    jbfill(bins,upper,lower,color(k,:),edge,0,0.5);
    hold on
    plot(bins,mn_rec, 'Color', color(k,:));

%     se_nrec=std(spk_dens_nrec_norm,0,1)./sqrt(numel(notrec_ids));
%     mn_nrec=mean(spk_dens_nrec_norm,1)+100;
%     upper=mn_nrec+se_nrec;
%     lower=mn_nrec-se_nrec;
%     color=[0 0 0];
%     edge=[0 0 0];
%     jbfill(bins,upper,lower,color,edge,0,0.5);
%     plot(bins,mn_nrec);
    title('Spike Density of recruited and not recruited neurons');
    xlim([1 2000]);
    ylim([-20 500]);
    % select recruited neurons
    %spikeTs_rec = spikeTs()
end
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


